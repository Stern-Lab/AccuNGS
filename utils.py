import os
import gzip
import shutil
import numpy as np

from Bio import SeqIO


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
