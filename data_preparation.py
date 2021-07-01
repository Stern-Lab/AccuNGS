"""
This is usually called by runner.py

Input: directory containing fastq/gz files or a directory containing such directories.
Output: fastq files in sizes ready for efficient processing.

Depending on currently availble RAM and CPUs or given values in -mm / --max_memory and -cc / --cpu_count the script
will divide the files into an efficient number of files to later run the processing script on.
If given an merge_opposing it will also merge the forward and backward reads in corresponding fastq files.
"""

import os
import gzip
import argparse
from functools import partial
import multiprocessing as mp

import psutil
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

from logger import pipeline_logger
from utils import get_files_by_extension, extract_gz


def merger_generator(forward_handle, reverse_handle, rep_length, log):
    for a, b in zip(SeqIO.parse(forward_handle, "fastq"), SeqIO.parse(reverse_handle, "fastq")):
        if a.id.split(" ")[0] != b.id.split(" ")[0]:
            # TODO: what does this mean and shouldnt it be an exception?
            log.warning("Problem, discrepancy in pair id's: {}, {}".format(a.id.split(" ")[0], b.id.split(" ")[0]))
        new_seq_id = a.id.split(" ")[0]
        new_seq_str = str(a.seq) + ("N" * rep_length) + str(b.seq)
        a_quals = a.letter_annotations["phred_quality"]
        b_quals = b.letter_annotations["phred_quality"]
        new_seq_qual = a_quals + [1.0 for a in range(rep_length)] + b_quals
        new_seq = SeqRecord(Seq.Seq(new_seq_str), id=new_seq_id, description="",
                            letter_annotations={"phred_quality": new_seq_qual})
        yield new_seq


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    Adjusted from  https://biopython.org/wiki/Split_large_file
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def get_fastq_records_num(fastq_file):
    i = 0
    with open(fastq_file) as f:
        for i, l in enumerate(f, 1):
            pass
    return i / 4


def split_fastq_file(fastq_file, output_dir, cpu_count, max_memory):
    fastq_records_num = get_fastq_records_num(fastq_file)
    part_size = fastq_records_num / cpu_count
    if not max_memory:
        max_memory = 0.9 * psutil.virtual_memory().available / 1000000  # 90% of currently available RAM.
    approx_memory_usage = part_size * cpu_count / 4  # rough estimate derived from testing
    while approx_memory_usage > max_memory:
        part_size /= 2
        approx_memory_usage = part_size * cpu_count / 4
    record_iter = SeqIO.parse(open(fastq_file), "fastq")
    fastq_file_name = os.path.basename(fastq_file)
    for i, batch in enumerate(batch_iterator(record_iter, part_size)):
        filename = os.path.join(output_dir, f"{fastq_file_name}.part_{i + 1}")
        with open(filename, "w") as handler:
            SeqIO.write(batch, handler, "fastq")


def find_read_files(input_dir, log):
    fastq_files = get_files_by_extension(input_dir, 'fastq')
    if len(fastq_files) == 0:
        gz_files = get_files_by_extension(input_dir, 'gz')
        file_type = "gz"
        files = gz_files
    else:
        log.debug(f"Found fastq files so not looking for gz files.")
        file_type = "fastq"
        files = fastq_files
    return files, file_type


def merge_opposing_reads(file1, file2, output_file, rep_length, file_type, log):
    read = partial(gzip.open, mode='rt') if file_type == 'gz' else partial(open, mode="r")
    write = partial(gzip.open, mode='wt') if file_type == 'gz' else partial(open, mode="w")
    if file_type == 'gz':
        output_file = output_file + ".gz"
    with read(file1) as forward_handle:
        with read(file2) as reverse_handle:
            with write(output_file) as merged_handle:
                SeqIO.write(merger_generator(forward_handle, reverse_handle, rep_length, log), merged_handle, "fastq")
    return output_file


def prepare_data_in_dir(input_dir, output_dir, rep_length, overlapping_reads, log, cpu_count, max_memory):
    input_dir_name = os.path.basename(input_dir)
    if input_dir_name == "":
        input_dir_name = os.path.basename(input_dir[:-1])  # that last '/' confuses basename..
    files, file_type = find_read_files(input_dir, log)
    if len(files) == 0:
        log.warning(f"Did not find any relevant files in {input_dir} !")
        return None
    if overlapping_reads == "Y" or overlapping_reads == "M":
        if len(files) == 2:
            log.debug(f"Found 2 opposing read files in {input_dir}")
            merged_reads = os.path.join(output_dir, input_dir_name + '_merged_reads.fastq')
            merged_reads = merge_opposing_reads(file1=files[0], file2=files[1], output_file=merged_reads,
                                                rep_length=rep_length, file_type=file_type, log=log)
            if file_type == 'gz':
                merged_reads = extract_gz(merged_reads, output_dir=output_dir)
            split_fastq_file(fastq_file=merged_reads, output_dir=output_dir, cpu_count=cpu_count, max_memory=max_memory)
        else:
            raise Exception(f"There are not 2 files to merge - overlapping or partial overlapping reads.")
    else:
        for file in files:
            if file_type == 'gz':
                file = extract_gz(file, output_dir=output_dir)
            split_fastq_file(fastq_file=file, output_dir=output_dir, cpu_count=cpu_count, max_memory=max_memory)


def prepare_data(input_dir, output_dir, cpu_count, max_memory, overlapping_reads, rep_length=60):
    log = pipeline_logger(logger_name='Data-Preparation', log_folder=output_dir)
    if not cpu_count:
        cpu_count = mp.cpu_count()
    os.makedirs(output_dir, exist_ok=True)
    prepare_data_in_dir(input_dir=input_dir, output_dir=output_dir, rep_length=rep_length, max_memory=max_memory,
                        overlapping_reads=overlapping_reads, log=log, cpu_count=cpu_count)
    sub_dirs = [f.path for f in os.scandir(input_dir) if f.is_dir()]
    for dir_path in sub_dirs:
        prepare_data_in_dir(input_dir=dir_path, output_dir=output_dir, rep_length=rep_length,
                            max_memory=max_memory, overlapping_reads=overlapping_reads, log=log, cpu_count=cpu_count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing fastq or gz files of one specific sample")
    parser.add_argument("-o", "--output_dir", help="Where the output files go", required=True)
    parser.add_argument("-r", "--rep_length", help="amount of N bases to repeat (default: 60)", type=int)
    parser.add_argument("-on", "--merge_opposing", default=None, nargs="+", type=str,
                        help="Notation of overlapping reads in the same directory to merge. Passing N would run without"
                             " considering overlapping reads (default is: '_R1 _R2')")
    parser.add_argument("-cc", "--cpu_count", default=None, type=int,
                        help="How many cpu's you have will determine to how many parts to split the file")
    args = parser.parse_args()
    prepare_data(input_dir=args.input_dir, output_dir=args.output_dir, rep_length=args.rep_length,
                 overlapping_reads=args.overlap_notation, cpu_count=args.cpu_count)
