import configparser
import os
import gzip
import shutil
import sys
import decimal
import hashlib
from _hashlib import HASH as Hash
from pathlib import Path
from typing import Union
from shutil import copyfile
import numpy as np
import pandas as pd

from Bio import SeqIO, Seq


def get_files_by_extension(dir_path, extension):
    files = []
    for filename in os.listdir(dir_path):
        if filename.endswith(f".{extension}"):
            files.append(os.path.join(dir_path, filename))
    return files


def get_files_in_dir(dir_path):
    files = []
    for filename in os.listdir(dir_path):
        files.append(os.path.join(dir_path, filename))
    return files


def extract_gz(gz_file, output_dir):
    file_name = os.path.basename(gz_file)
    output_file = os.path.join(output_dir, file_name[:-3])
    with gzip.open(gz_file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return output_file


def is_any_nan(var):
    return bool(var is None) or bool(np.isnan(var))


def drange(x, y, jump):
    """Adapted from https://stackoverflow.com/questions/7267226/range-for-floats"""
    jump = str(jump)  # necessary to avoid floating point chaos
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)


def concatenate_files_by_extension(input_dir, extension, output_path, remove_headers=True):
    """
    Concatenates all files in a directory without loading them all into memory.
    By default skip first line of all files except the first one to avoid multiple headers.
    """
    files = get_files_by_extension(input_dir, extension)
    is_first_file = True
    with open(output_path, "w") as output_handle:
        for input_file in files:
            with open(input_file, "r") as input_handle:
                is_first_line = True
                for line in input_handle:
                    if is_first_line and not is_first_file and remove_headers:
                        # skip the first line
                        is_first_line = False
                    else:
                        output_handle.write(line)
                is_first_file = False


def get_sequence_from_fasta(fasta_file):
    with open(fasta_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) != 1:
        raise Exception("fasta file must contain only one record!")
    return str(records[0].seq)


def md5_update_from_file(filename: Union[str, Path], hash: Hash) -> Hash:
    assert Path(filename).is_file(), f"{filename} is not a valid file!"
    with open(str(filename), "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash.update(chunk)
    return hash


def md5_file(filename: Union[str, Path]) -> str:
    return str(md5_update_from_file(filename, hashlib.md5()).hexdigest())


def md5_update_from_dir(directory: Union[str, Path], hash: Hash) -> Hash:
    assert Path(directory).is_dir(), f"{directory} is not a valid directory!"
    for path in sorted(Path(directory).iterdir(), key=lambda p: str(p).lower()):
        hash.update(path.name.encode())
        if path.is_file():
            hash = md5_update_from_file(path, hash)
        elif path.is_dir():
            hash = md5_update_from_dir(path, hash)
    return hash


def md5_dir(directory: Union[str, Path]) -> str:
    # https://stackoverflow.com/questions/24937495/how-can-i-calculate-a-hash-for-a-filesystem-directory-using-python
    return str(md5_update_from_dir(directory, hashlib.md5()).hexdigest())


def create_new_ref_with_freqs(reference_fasta_file, freqs_file, min_coverage, output_file, drop_indels, freq_threshold):
    # TODO: what about deletions in the start or begining?
    """Create reference from freqs filling unaligned parts with the given reference file."""
    with open(reference_fasta_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) != 1:
        raise Exception("fasta file must contain only one record!")
    record = records[0]

    freqs_df = pd.read_table(freqs_file)
    freqs_df = freqs_df[freqs_df["base_rank"] == 0]
    freqs_df.loc[freqs_df["coverage"] <= min_coverage, 'read_base'] = np.nan
    freqs_df = freqs_df[freqs_df['frequency'] > float(freq_threshold)]

    if drop_indels:
        freqs_df = freqs_df[freqs_df["ref_pos"] == np.round(freqs_df['ref_pos'])]  # drop insertions
        freqs_df = freqs_df[freqs_df["read_base"] != "-"]  # drop deletions

    # create the new conc - not depending on the reference
    new_seq = freqs_df['read_base'].fillna('N')
    new_seq = new_seq.replace('-', np.nan).dropna()  # drop deletions
    new_seq = new_seq.str.cat()
    new_sequence = Seq.Seq(new_seq)
    record.seq = new_sequence
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")


def get_mp_results_and_report(async_objects_list):
    """Expects a list of pool.apply_async objects, runs them and reports when each of them is done."""
    i = 1
    results = []
    for async_object in async_objects_list:
        try:
            res = async_object.get()
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise
        results.append(res)
        sys.stdout.write(f"\rDone {i}/{len(async_objects_list)} parts.")
        sys.stdout.flush()
        i = i + 1
    sys.stdout.write("\n")
    return results


def reverse_string(string):
    # From here: https://stackoverflow.com/questions/931092/reverse-a-string-in-python
    return string[::-1]


def get_config():
    accungs_dir = os.path.dirname(os.path.abspath(__file__))
    config_file = os.path.join(accungs_dir, 'config.ini')
    if not os.path.isfile(config_file):
        sample_config = os.path.join(accungs_dir, "config_sample.ini")
        copyfile(sample_config, config_file)
    config = configparser.ConfigParser()
    config.read(config_file)
    return config


def create_fixed_param_list(input_dir, output_dir, params_dict):
    """
    Use to created a parameter list for pbs_multi_runner when you have different inputs in sub directories of
    input_dir dirs and want to run the same parameters on all of them.
    """
    input_dirs = [d for d in get_files_in_dir(input_dir) if os.path.isdir(d)]
    dirs_dict = {in_dir: os.path.join(output_dir, os.path.basename(in_dir)) for in_dir in input_dirs}
    param_list = []
    for in_dir, out_dir in dirs_dict.items():
        param = {'input_dir': in_dir, 'output_dir': out_dir}
        param.update(params_dict)
        param_list.append(param)
    return param_list
