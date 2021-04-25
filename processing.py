"""
This file is used by runner.py for base calling on each of the parts of the input files created by
data_preparation.py (parallelization is done in runner.py).

Input: fastq files from part I and a reference fasta file.
Output: called_bases - every base called and its attributes.
        ignored_bases - every base ignored, its attributes and why it was ignored.
        suspicious_reads - every suspicious read and why it is suspicious.
        ignored_reads - every read ignored and why.
        read_counter - every read and a number describing how many times it was aligned by blast.
        blast files - alignment output files.

This is where the main logic of the pipeline happens and the runner runs this on each of the fastq files in parallel.
Basically it first aligns the reads to the reference with blast and then goes over every read and nucleotide and decides
whether to filter them out or leave them in.
"""

import argparse
import os

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from logger import pipeline_logger
import numpy as np
from utils import reverse_string


def run_blast(reads_fasta, reference, output, mode, task, evalue, perc_identity, num_alignments, dust, soft_masking,
              log):
    if mode == "SeqToRef":
        query = reads_fasta
        subject = reference
    elif mode == "RefToSeq":
        make_db = NcbimakeblastdbCommandline(input_file=reads_fasta, dbtype='nucl')
        stdout, stderr = make_db()
        for output_string in [stdout, stderr]:
            if output_string:
                log.debug(output_string)
        query = reference
        subject = reads_fasta
    else:
        raise ValueError(f"parameter mode must be one of ['RefToSeq', 'SeqToRef'] and not '{mode}' ! ")
    outfmt = "6 qseqid sseqid qstart qend qstrand sstart send sstrand length btop qseq sseq"  # defines blast output
    blast_instance = NcbiblastnCommandline(query=query, subject=subject, task=task, out=output, dust=dust,
                                           num_alignments=num_alignments, soft_masking=soft_masking,
                                           perc_identity=perc_identity, evalue=evalue, outfmt=outfmt)
    stdout, stderr = blast_instance()
    return stdout, stderr


def _rename_columns(df, mode):
    if mode == 'SeqToRef':
        new_name = {'query': 'read', 'subject': 'ref'}
    else:
        new_name = {'query': 'ref', 'subject': 'read'}
    new_columns = {name: name.replace('query', new_name['query']).replace('subject', new_name['subject']) for name in
                   df.columns}
    return df.rename(columns=new_columns)


def get_alignments(blast_output, fastq_file, reads_overlap, mode):
    df = pd.read_csv(blast_output, sep="\t", index_col=False,
                     names=["query_id", "subject_id", "query_start", 'query_end', 'subject_start',
                            'subject_end', 'plus_or_minus', 'length', 'mutations', 'query_seq', 'subject_seq'])
    alignments = _rename_columns(df, mode)
    suspicious_reads = pd.DataFrame()
    alignments, ignored_reads, multi_mapped_alignments, read_counter = filter_reads_by_alignment_count(alignments,
                                                                                                       reads_overlap)
    multi_mapped_alignments['suspicious_because'] = "multiple alignments"
    suspicious_reads = suspicious_reads.append(multi_mapped_alignments)
    quality = get_quality(fastq_file, alignments)
    alignments['quality'] = alignments.read_id.map(lambda r: quality[r])
    return alignments, ignored_reads.reset_index(drop=True), suspicious_reads.reset_index(drop=True), read_counter


def filter_reads_by_alignment_count(alignments, reads_overlap):
    ignored_reads = pd.DataFrame()
    read_counter = alignments.groupby('read_id').size()
    multiple_reads = read_counter[read_counter > 2].index
    multiple_alignments = alignments[alignments['read_id'].isin(multiple_reads)].copy()
    if reads_overlap:
        single_reads = list(read_counter[read_counter == 1].index)
        ignored_reads['read_id'] = single_reads
        ignored_reads['dropped_because'] = "single read"
        alignments = alignments.query("read_id not in @single_reads")
    read_counter.name = "number_of_alignments"
    return alignments, ignored_reads, multiple_alignments, read_counter


def get_quality(fastq_file, alignments):
    quality = {}
    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in alignments.read_id.values:
                quality[record.id] = record.letter_annotations["phred_quality"]
    return quality


def get_max_insertion_value(df_index, insertion):
    max_decimal = 0
    for num in range(10):
        decimal = num / 10
        if insertion + decimal in df_index:
            if decimal > max_decimal:
                max_decimal = decimal
    return insertion + max_decimal


def fix_insertions_index(df, ref_start):
    # change insertions index to original index + 0.001
    df.index = df.index + ref_start
    insertions = df[df.ref_seq == "-"].index.to_list()
    df_index = [float(i) for i in df.index.to_list()]
    insertion_counter = 0
    for insertion in insertions:
        insertion = insertion - insertion_counter  # each insertion moves the other insertions
        insertion_position = df_index.index(insertion)
        insertion_value = round(insertion + 0.001, 3)  # fix floating point nonsense
        while insertion_value - 1 in df_index:
            insertion_value = round(insertion_value + 0.001, 3)
        df_index[insertion_position] = insertion_value
        df_index = [x - 1 if x > insertion else x for x in df_index]
        insertion_counter += 1
    df.index = df_index
    return df


def parse_alignment_data(data, read_seq, ref_seq, read_start, ref_start, mode, minus_seq=False):
    if minus_seq and (mode == "SeqToRef"):
        ref_start = ref_start.replace('start', 'end')
        complimentary_base = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}
        seq = {'read': read_seq, 'ref': ref_seq}
        for seq_type, seq_data in seq.items():  # reverse and take complimentary
            seq_data = list(reverse_string(data[seq_data]))
            seq_data = [complimentary_base[base] for base in seq_data]
            seq[seq_type] = seq_data
    else:
        seq = {'read': list(data[read_seq]), 'ref': list(data[ref_seq])}
    if minus_seq and not (mode == "SeqToRef"):
        read_start = read_start.replace('start', 'end')
    df = pd.DataFrame(data={'read_seq': seq['read'],
                            'ref_seq': seq['ref']})
    df['read_pos'] = (df.index + data[read_start]).astype(int)
    quality = data['quality'].copy()
    deletions = df[df.read_seq == "-"].read_pos
    if len(deletions) > 0:
        deletions = df[df.read_seq == "-"].read_pos
        for deletion in deletions.values:  # deletions have no quality...!
            quality.insert(deletion - 1, np.inf)  # set deletions quality to be inf so that we don't filter them later
    df['quality'] = df.read_pos.map(lambda pos: quality[pos - 1])
    df = fix_insertions_index(df, data[ref_start])
    return df


def get_alignment_df(alignment_data, mode):
    minus_seq = False
    if alignment_data['plus_or_minus'] == 'minus':
        minus_seq = True
    bases = parse_alignment_data(alignment_data, read_seq='read_seq', ref_seq='ref_seq', read_start='read_start',
                                 ref_start='ref_start', mode=mode, minus_seq=minus_seq)
    bases['read_id'] = alignment_data['read_id']
    bases['alignment_id'] = alignment_data.name
    bases['ref_pos'] = bases.index
    bases['plus_or_minus'] = alignment_data['plus_or_minus']
    bases = bases.rename(columns={'read_seq': 'read_base', 'ref_seq': 'ref_base'})
    return bases.reset_index(drop=True)


def filter_target_mean_by(df, target_column, by, min_mean):
    target_nunique_values = df.groupby(by)[target_column].mean()
    target_nunique_values.name = 'tmpCol'
    target_nunique_values = target_nunique_values.reset_index()
    df = df.merge(target_nunique_values, on=by)
    df_with_high_mean = df[df['tmpCol'] >= min_mean].drop(columns=['tmpCol']).copy()
    df_with_low_mean = df[df['tmpCol'] < min_mean].drop(columns=['tmpCol']).copy()
    return df_with_high_mean, df_with_low_mean


def filter_target_nunique_by(df, target_column, by, required_value=1):
    target_nunique_values = df.groupby(by)[target_column].nunique()
    target_nunique_values.name = 'tmpCol'
    target_nunique_values = target_nunique_values.reset_index()
    df = df.merge(target_nunique_values, on=by)
    df_with_right_number_of_targets = df[df['tmpCol'] == required_value].drop(columns=['tmpCol']).copy()
    df_with_wrong_number_of_targets = df[df['tmpCol'] != required_value].drop(columns=['tmpCol']).copy()
    return df_with_right_number_of_targets, df_with_wrong_number_of_targets


def filter_bases(called_bases, quality_threshold, reads_overlap):
    called_bases, multi_aligned_bases = filter_target_nunique_by(called_bases, by=['read_id', 'read_pos'],
                                                                 target_column='ref_pos')
    multi_aligned_bases['dropped_because'] = "read position aligned to more than one ref position"
    called_bases, mismatching_bases = filter_target_nunique_by(called_bases, by=['read_id', 'ref_pos'],
                                                               target_column='read_base')
    mismatching_bases['dropped_because'] = "overlapping reads mismatch"
    called_bases, low_quality_bases = filter_target_mean_by(called_bases, by=['read_id', 'ref_pos'],
                                                            target_column='quality', min_mean=quality_threshold)
    low_quality_bases['dropped_because'] = f"base phred score lower than threshold: {quality_threshold}"
    ignored_bases = [multi_aligned_bases, mismatching_bases, low_quality_bases]
    if reads_overlap:
        called_bases, no_plus_minus_bases = filter_target_nunique_by(called_bases, by=['read_id', 'ref_pos'],
                                                                     target_column='plus_or_minus', required_value=2)
        no_plus_minus_bases['dropped_because'] = "did not align with both plus and minus"
        ignored_bases.append(no_plus_minus_bases)
    ignored_bases = pd.concat(ignored_bases)
    return called_bases, ignored_bases


def basecall(blast_output_file, fastq_file, output_dir, quality_threshold, mode, reads_overlap):
    base_filename = os.path.basename(fastq_file)
    alignments, ignored_reads, suspicious_reads, read_counter = get_alignments(blast_output_file, fastq_file=fastq_file,
                                                                               reads_overlap=reads_overlap, mode=mode)
    read_counter.to_csv(os.path.join(output_dir, base_filename + ".read_counter"), sep="\t")
    suspicious_reads.to_csv(os.path.join(output_dir, base_filename + ".suspicious_reads"), sep="\t")
    ignored_reads.to_csv(os.path.join(output_dir, base_filename + ".ignored_reads"), sep="\t")
    called_bases = alignments.apply(lambda row: get_alignment_df(row, mode=mode), axis=1)
    called_bases = pd.concat(list(called_bases)).reset_index(drop=True)
    called_bases, ignored_bases = filter_bases(called_bases, quality_threshold, reads_overlap)
    ignored_bases.to_csv(os.path.join(output_dir, base_filename + ".ignored_bases"), sep="\t", index=False)
    called_bases.to_csv(os.path.join(output_dir, base_filename + ".called_bases"), sep="\t", index=False)


def convert_fastq_to_fasta(output_dir, fastq_file):
    reads_fasta_file_name = os.path.basename(fastq_file).replace('fastq', 'fasta')
    reads_fasta_file_path = os.path.join(output_dir, reads_fasta_file_name)
    with open(fastq_file, "r") as input_handle:
        with open(reads_fasta_file_path, "w") as output_handle:
            SeqIO.convert(input_handle, "fastq", output_handle, "fasta")
    return reads_fasta_file_path


def process_fastq(fastq_file, reference, output_dir, quality_threshold, task, evalue, dust, num_alignments,
                  soft_masking, perc_identity, mode, reads_overlap):
    log = pipeline_logger(logger_name=f"Computation_{os.path.basename(fastq_file)}", log_folder=output_dir)
    blast_output = os.path.join(output_dir, 'blast')
    os.makedirs(blast_output, exist_ok=True)
    reads_fasta_file_path = convert_fastq_to_fasta(fastq_file=fastq_file, output_dir=blast_output)
    blast_output_file = reads_fasta_file_path + ".blast"
    run_blast(reads_fasta=reads_fasta_file_path, reference=reference, output=blast_output_file, mode=mode, task=task,
              evalue=evalue, num_alignments=num_alignments, dust=dust, soft_masking=soft_masking, log=log,
              perc_identity=perc_identity)
    basecall_output = os.path.join(output_dir, 'basecall')
    if not quality_threshold:
        quality_threshold = 30
    os.makedirs(basecall_output, exist_ok=True)
    basecall(blast_output_file=blast_output_file, fastq_file=fastq_file, output_dir=basecall_output,
             quality_threshold=quality_threshold, mode=mode, reads_overlap=reads_overlap)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq_file", required=True,
                        help="Path to fastq file containing reads")
    parser.add_argument("-r", "--reference_path", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-ro", "--reads_overlap", required=True)
    parser.add_argument("-q", "--quality_threshold", type=int, help="bases with quality under this will be ignored")
    parser.add_argument("-bt", "--blast_task", help="blast's task parameter")
    parser.add_argument("-be", "--blast_evalue", help="blast's evalue parameter", type=float)
    parser.add_argument("-bd", "--blast_dust", help="blast's dust parameter")
    parser.add_argument("-bn", "--blast_num_alignments", type=int, help="blast's num_alignments parameter")
    parser.add_argument("-bp", "--blast_perc_identity", type=int, help="blast's perc_identity parameter")
    parser.add_argument("-bs", "--blast_soft_masking", type=int, help="blast's soft_masking parameter")
    parser.add_argument("-bm", "--blast_mode", help="RefToSeq or SeqToRef")
    args = parser.parse_args()
    process_fastq(fastq_file=args.fastq_file, reference=args.reference_path, output_dir=args.output_dir,
                  reads_overlap=args.reads_overlap, quality_threshold=args.quality_threshold, task=args.blast_task,
                  evalue=args.blast_evalue, dust=args.blast_dust, num_alignments=args.blast_num_alignments,
                  mode=args.blast_mode, perc_identity=args.blast_perc_identity, soft_masking=args.blast_soft_masking)
