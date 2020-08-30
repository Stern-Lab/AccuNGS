import os
import gzip
import shutil
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


def extract_gz(gz_file):
    output_file = gz_file[:-3]
    with gzip.open(gz_file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return output_file


def is_any_nan(var):
    return bool(var is None) or bool(np.isnan(var))


def concatenate_files_by_extension(input_dir, extension, output_path, remove_duplicate_headers=True):
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
                    if is_first_line and not is_first_file and remove_duplicate_headers:
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


def create_new_ref_with_freqs(reference_fasta_file, freqs_file, min_coverage, output_file, drop_indels):
    # TODO: what about deletions in the start or begining?
    """Create reference from freqs filling unaligned parts with the given reference file."""
    with open(reference_fasta_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) != 1:
        raise Exception("fasta file must contain only one record!")
    record = records[0]
    ref = pd.DataFrame(list(str(record.seq)), columns=['ref_base_from_fasta'])
    ref.index = (ref.index + 1).astype(float)
    df = pd.read_table(freqs_file)
    df.loc[df["coverage"] >= min_coverage, 'coverage'] = None
    df = df[df["base_rank"] == 0]
    if drop_indels:
        df = df[df["ref_pos"] == np.round(df['ref_pos'])]   # drop insertions
        df = df[df["read_base"] != "-"]                     # drop deletions
    df = df.merge(ref, left_on='ref_pos', right_index=True, how='outer').sort_values(by='ref_pos')
    new_seq = df['read_base'].fillna(df['ref_base_from_fasta']).replace('-', np.nan).dropna().str.cat()
    new_sequence = Seq.Seq(new_seq)
    record.seq = new_sequence
    with open(output_file, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
