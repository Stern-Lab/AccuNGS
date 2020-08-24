"""
This script accepts 1 or 2 fastq / gz files, concatenates them and extracts them if necessary and then splits them into
several smaller files ready for parallel execution.
"""

import os
import gzip
import argparse
from functools import partial
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

from logger import pipeline_logger
from utils import get_files_by_extension, extract_gz


def merger_generator(forward_handle,reverse_handle, rep_length, log):
    # TODO: why not replace mismatching bases with N?
    for a, b in zip(SeqIO.parse(forward_handle, "fastq"), SeqIO.parse(reverse_handle, "fastq")):
        if a.id.split(" ")[0] != b.id.split(" ")[0]:
            # TODO: what does this mean and shouldnt it be an exception?
            log.warning("Problem, discrepancy in pair id's: {}, {}".format(a.id.split(" ")[0], b.id.split(" ")[0]))
        new_seq_id = a.id.split(" ")[0]
        new_seq_str = str(a.seq) + ("N"*rep_length) + str(b.seq)
        a_quals = a.letter_annotations["phred_quality"]
        b_quals = b.letter_annotations["phred_quality"]
        new_seq_qual = a_quals+[1.0 for a in range(rep_length)]+b_quals
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


def split_fastq_file(fastq_file, output_dir, part_size):
    record_iter = SeqIO.parse(open(fastq_file), "fastq")
    fastq_file_name = os.path.basename(fastq_file)
    for i, batch in enumerate(batch_iterator(record_iter, part_size)):  # TODO: parallelize this?
        filename = os.path.join(output_dir, f"{fastq_file_name}.part_{i+1}")
        with open(filename, "w") as handler:
            SeqIO.write(batch, handler, "fastq")


def find_relevant_files(input_dir, log):
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


def are_opposing(files, opposing_strings=None):
    if (opposing_strings[0] in files[0] and opposing_strings[1] in files[1]) or (
            opposing_strings[1] in files[0] and opposing_strings[0] in files[1]):
        return True
    else:
        return False


def prepare_data_in_dir(input_dir, output_dir, part_size, rep_length, opposing_strings, log):
    input_dir_name = os.path.basename(input_dir)
    if input_dir_name == "":
        input_dir_name = os.path.basename(input_dir[:-1])  # that last '/' confuses basename..
    files, file_type = find_relevant_files(input_dir, log)
    if len(files) == 0:
        log.warning(f"Did not find any relevant files in {input_dir} !")
        return None
    if are_opposing(files, opposing_strings):
        if len(files) == 2:
            log.debug(f"Found 2 opposing files in {input_dir}")
            merged_reads = os.path.join(output_dir, input_dir_name + '_merged_reads.fastq')
            merged_reads = merge_opposing_reads(file1=files[0], file2=files[1], output_file=merged_reads,
                                                rep_length=rep_length, file_type=file_type, log=log)
            if file_type == 'gz':
                merged_reads = extract_gz(merged_reads)
            split_fastq_file(fastq_file=merged_reads, output_dir=output_dir, part_size=part_size)
        else:
            raise Exception(f"Found more than 2 files containing opposing_strings: {opposing_strings} !")
    else:
        for file in files:
            if file_type == 'gz':
                file = extract_gz(file)
            split_fastq_file(fastq_file=file, output_dir=output_dir, part_size=part_size)


def prepare_data(input_dir, output_dir, rep_length=None, part_size=None, opposing_strings=None):
    log = pipeline_logger(logger_name='Data-Preparation', log_folder=output_dir)
    if rep_length is None:
        rep_length = 60
    if part_size is None:
        part_size = 1000
    if opposing_strings is None:
        opposing_strings = ("_R1", "_R2")
    os.makedirs(output_dir, exist_ok=True)
    prepare_data_in_dir(input_dir=input_dir, output_dir=output_dir, part_size=part_size, rep_length=rep_length,
                        opposing_strings=opposing_strings, log=log)
    sub_dirs = [f.path for f in os.scandir(input_dir) if f.is_dir()]
    for dir_path in sub_dirs:
        prepare_data_in_dir(input_dir=dir_path, output_dir=output_dir, part_size=part_size, rep_length=rep_length,
                            opposing_strings=opposing_strings, log=log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing fastq or gz files of one specific sample")
    parser.add_argument("-o", "--output_dir", help="Where the output files go", required=True)
    parser.add_argument("-p", "--part_size", type=int,
                        help="Number of parts to split fastq file into as preparation for parallel execution. "
                             "(default: 1000)")
    parser.add_argument("-r", "--rep_length", help="amount of N bases to repeat (default: 60)", type=int) #TODO: is this the right default?
    parser.add_argument("-op", "--opposing_strings", default=None, type=str,
                        help="A tuple with strings which represent the forward and backward reads "
                             "(default is: ('R1','R2'))")
    args = parser.parse_args()
    prepare_data(input_dir=args.input_dir, output_dir=args.output_dir, rep_length=args.rep_length,
                 part_size=args.part_size, opposing_strings=args.opposing_strings)
