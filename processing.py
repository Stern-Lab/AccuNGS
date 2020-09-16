import argparse
import os
from decimal import Decimal

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from logger import pipeline_logger
import numpy as np


def run_blast(reads_fasta, reference, output, mode, task, evalue, perc_identity, num_alignments, dust, soft_masking, log):
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
    outfmt = '6 qseqid sseqid qstart qend qstrand sstart send sstrand length btop qseq sseq'  # defines blast output
    blast_instance = NcbiblastnCommandline(query=query, subject=subject, task=task, out=output, dust=dust,
                                           num_alignments=num_alignments, soft_masking=soft_masking,
                                           perc_identity=perc_identity, evalue=evalue, outfmt=outfmt)
    stdout, stderr = blast_instance()
    return stdout, stderr


def _rename_columns(df):
    new_name = {}
    for col in ['subject_id', 'query_id']:
        if df[col].nunique() == 1:
            new_name[col[:-3]] = 'ref'
        else:
            new_name[col[:-3]] = 'read'
    if new_name['subject'] == new_name['query']:
        raise Exception("Blast returned 2 constant columns! This really shouldn't happen.")
    new_columns = {name: name.replace('query', new_name['query']).replace('subject', new_name['subject']) for name in df.columns}
    return df.rename(columns=new_columns)


def get_alignments(blast_output, fastq_file, reads_overlap):
    df = pd.read_csv(blast_output, sep="\t", index_col=False,
                     names=["query_id", "subject_id", "query_start", 'query_end', 'subject_start',
                            'subject_end', 'plus_or_minus', 'length', 'mutations',  'query_seq', 'subject_seq'])
    alignments = _rename_columns(df)
    suspicious_reads = pd.DataFrame()
    alignments, ignored_reads, multi_mapped_alignments, read_counter = filter_reads_by_alignment_count(alignments,
                                                                                                       reads_overlap)
    multi_mapped_alignments['suspicious_because'] = "multiple alignments"
    suspicious_reads = suspicious_reads.append(multi_mapped_alignments)
    if reads_overlap:
        alignments, multi_paired_alignments = align_pairs(alignments)
        multi_paired_alignments['suspicious_because'] = "aligned with more than one pair"
        suspicious_reads = suspicious_reads.append(multi_paired_alignments)
    quality = get_quality(fastq_file, alignments)
    alignments['quality'] = alignments.read_id.map(lambda r: quality[r])

    #TODO: drop multi aligned bases without overlap
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


def _apply_find_minus_with_max_range(row, minuses):
    max_intersection_len = 0
    max_i = None
    plus_range = range(row.ref_start, row.ref_end+1)
    for i, minus_row in minuses.iterrows():
        if minus_row.ref_end < minus_row.ref_start:
            minus_range = range(minus_row.ref_end, minus_row.ref_start+1)
        else:
            minus_range = range(minus_row.ref_start, minus_row.ref_end+1)
        intersection = set(minus_range) & set(plus_range)
        if len(intersection) > max_intersection_len:
            max_intersection_len = len(intersection)
            max_i = i
    if max_i:
        return {'pair_index': max_i, 'intersection_length': max_intersection_len}
    else:
        return {'pair_index': None, 'intersection_length': None}


def _find_pairs(df):
    minuses = df[df.plus_or_minus == 'minus']
    pluses = df[df.plus_or_minus == 'plus']
    if pluses.empty:
        return None
    df_with_intersection = pluses.apply(lambda row: _apply_find_minus_with_max_range(row, minuses), axis=1,
                                        result_type='expand')
    return df_with_intersection


def get_minus_values(data, values):
    df = data.copy()
    values = values + ['pair']
    df = df.merge(df[values], left_index=True, right_on='pair', suffixes=['_minus', '_plus'], how='inner')
    return df


def get_quality(fastq_file, alignments):
    quality = {}
    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in alignments.read_id.values:
                quality[record.id] = record.letter_annotations["phred_quality"]
    return quality


def get_paired_info(df):
    # find pairs of plus and minus of the same read and return the longest alignment
    paired_info = df.groupby('read_id', sort=False).apply(lambda read_df: _find_pairs(read_df))
    # return to df's index
    paired_info = paired_info.reset_index().set_index('level_1')
    return paired_info['pair_index']


def align_pairs(df):
    df = df.copy()
    paired_info = get_paired_info(df)
    df['pair'] = paired_info
    df = get_minus_values(df, ['ref_seq', 'read_seq', 'ref_start', 'ref_end', 'read_start', 'read_end'])
    multi_paired_alignments = df.dropna(subset=['pair'])
    multi_paired_alignments = multi_paired_alignments[multi_paired_alignments['pair'].duplicated()]
    just_minus = df[df.plus_or_minus == 'minus'].copy()  # We put all the relevant data in the minus reads
    alignments = just_minus.reset_index(drop=True)
    alignments['multi_mapped_bases'] = get_multi_mapped_bases_wrapper(alignments)
    return alignments, multi_paired_alignments


def _get_overlapping_bases_on_other_rows(this_row, start_col, end_col, df):
    this_range = range(this_row[start_col], this_row[end_col])
    total_intersections = set()
    for i, that_row in df.iterrows():
        that_range = range(that_row[start_col], that_row[end_col])
        intersection = set(this_range) & set(that_range)
        total_intersections = total_intersections | intersection
    if len(total_intersections) == 0:
        return set()
    return total_intersections


def _get_multimapped_bases_on_plus_and_minus(this_row, df):
    plus_multimapped_bases = _get_multimapped_bases(this_row, start_col='read_start_plus', end_col='read_end_plus',
                                                    df=df)
    minus_multimapped_bases = _get_multimapped_bases(this_row, start_col='read_start_minus', end_col='read_end_minus',
                                                     df=df)
    return plus_multimapped_bases | minus_multimapped_bases


def _get_multimapped_bases(this_row, start_col, end_col, df):
    this_range = range(this_row[start_col], this_row[end_col]+1)
    total_intersections = set()
    for i, that_row in df.iterrows():
        that_range = range(that_row[start_col], that_row[end_col]+1)
        intersection = set(this_range) & set(that_range)
        total_intersections = total_intersections | intersection
    if len(total_intersections) == 0:
        return set()
    return total_intersections


def get_multi_mapped_bases_wrapper(data):
    """ignore bases which were aligned to more than one pair"""
    int_cols = data[['read_start_plus', 'read_end_plus', 'read_start_minus', 'read_end_minus']].astype(int)
    rel_cols = int_cols.join(data['read_id'])
    multi_mapped_bases_per_row = rel_cols.groupby('read_id').apply(
        lambda df: df.apply(
            lambda row: _get_multimapped_bases_on_plus_and_minus(row, df.drop(row.name)), axis=1))
    return multi_mapped_bases_per_row.reset_index().set_index('level_1')[0]


def reverse_string(string):
    # From here: https://stackoverflow.com/questions/931092/reverse-a-string-in-python
    return string[::-1]


def get_max_insertion_value(df_index, insertion):
    max_decimal = 0
    for num in range(10):
        decimal = num/10
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
        insertion = insertion - insertion_counter                   # each insertion moves the other insertions
        insertion_position = df_index.index(insertion)
        insertion_value = round(insertion + 0.001, 3)               # fix floating point nonsense
        while insertion_value-1 in df_index:
            insertion_value = round(insertion_value + 0.001, 3)
        df_index[insertion_position] = insertion_value
        df_index = [x - 1 if x > insertion else x for x in df_index]
        insertion_counter += 1
    df.index = df_index
    return df


def create_read_df(data, read_seq, ref_seq, read_start, ref_start, mode, minus_seq=False):
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
        for deletion in deletions.values:         # deletions have no quality...!
            quality.insert(deletion - 1, np.inf)  # set deletions quality to be inf so that we dont filter them later
    df['quality'] = df.read_pos.map(lambda pos: quality[pos - 1])
    df = fix_insertions_index(df, data[ref_start])
    return df


def create_plus_minus_dfs(data, mode):
    plus_df = create_read_df(data, read_seq='read_seq_plus', ref_seq='ref_seq_plus', read_start='read_start_plus',
                             ref_start='ref_start_plus', mode=mode)
    minus_df = create_read_df(data, read_seq='read_seq_minus', ref_seq='ref_seq_minus', read_start='read_start_minus',
                              ref_start='ref_start_minus', minus_seq=True, mode=mode)
    return plus_df, minus_df


def filter_mismatching_bases(df):
    mismatching_bases = df[df['read_seq_plus'] != df['read_seq_minus']].copy()
    mismatching_bases['dropped_because'] = "paired end bases mismatch"
    df = df[~df.index.isin(mismatching_bases.index)].copy()
    return df, mismatching_bases


def filter_overlapping_bases(df, multimapped_bases):
    multi_mapped_bases = df[df.read_pos_minus.isin(multimapped_bases) |
                            df.read_pos_plus.isin(multimapped_bases)].copy()
    multi_mapped_bases['dropped_because'] = "multimapped bases"
    df = df[~df.index.isin(multi_mapped_bases.index)].copy()
    return df, multi_mapped_bases


def filter_low_quality_bases(data, quality_threshold, quality_column):
    df = data.copy()
    low_quality_bases = df[df[quality_column] < quality_threshold].copy()
    low_quality_bases['dropped_because'] = f"average quality under quality threshold: {quality_threshold}"
    df = df[~df.index.isin(low_quality_bases.index)]
    return df, low_quality_bases


def basecall_on_overlapping_reads(read_data, quality_threshold, mode):
    plus_df, minus_df = create_plus_minus_dfs(read_data, mode)
    bases = plus_df.join(minus_df, how='outer', lsuffix='_plus', rsuffix='_minus')
    bases['read_id'] = read_data['read_id']
    bases['avg_quality'] = (bases['quality_plus'] + bases['quality_minus']) / 2
    bases, mismatching_bases = filter_mismatching_bases(bases)
    bases, low_quality_bases = filter_low_quality_bases(bases, quality_threshold, quality_column='avg_quality')
    dropped_bases = pd.concat([mismatching_bases, low_quality_bases])
    bases = bases.rename(columns={'read_seq_plus': 'read_base', 'ref_seq_plus': 'ref_base'}).drop(
        ['read_seq_minus', 'ref_seq_minus', 'avg_quality'], axis=1)  # these are redundant columns at this point.
    return bases, dropped_bases


def basecall_on_read_no_overlap(read_data, quality_threshold, mode):
    minus_seq = False
    if read_data['plus_or_minus'] == 'minus':
        minus_seq = True
    bases = create_read_df(read_data, read_seq='read_seq', ref_seq='ref_seq', read_start='read_start',
                           ref_start='ref_start', mode=mode, minus_seq=minus_seq)
    bases['read_id'] = read_data['read_id']
    bases, dropped_bases = filter_low_quality_bases(bases, quality_threshold, quality_column='quality')
    bases = bases.rename(columns={'read_seq': 'read_base', 'ref_seq': 'ref_base'})
    return bases, dropped_bases


def filter_multimapped_bases(called_bases, ignored_bases, reads_overlap):
    if reads_overlap:
        columns = ['read_pos_minus', 'read_pos_plus']
    else:
        columns = ['read_pos']
    multimapped_bases = []
    for col in columns:
        multimapped_bases.append(called_bases[called_bases.duplicated([col, 'read_id'], keep=False)])
    multimapped_bases = pd.concat(multimapped_bases)
    multimapped_bases['dropped_because'] = "multimapped read position"
    ignored_bases = pd.concat([ignored_bases, multimapped_bases])
    for col in columns:
        called_bases = called_bases.drop_duplicates([col, 'read_id'], keep=False)
    return called_bases, ignored_bases


def basecall(blast_output_file, fastq_file, output_dir, quality_threshold, mode, reads_overlap):
    """
    outputs ignored_reads, ignored_bases, suspicious_alignments, called_bases
    """
    base_filename = os.path.basename(fastq_file)
    alignments, ignored_reads, suspicious_reads, read_counter = get_alignments(blast_output_file, fastq_file=fastq_file, reads_overlap=reads_overlap)
    read_counter.to_csv(os.path.join(output_dir, base_filename + ".read_counter"), sep="\t")
    suspicious_reads.to_csv(os.path.join(output_dir, base_filename+".suspicious_reads"), sep="\t")
    ignored_reads.to_csv(os.path.join(output_dir, base_filename+".ignored_reads"), sep="\t")
    if reads_overlap:
        bases_object = alignments.apply(
            lambda row: basecall_on_overlapping_reads(row, quality_threshold=quality_threshold, mode=mode),
            axis=1, result_type="expand")
    else:
        bases_object = alignments.apply(
            lambda row: basecall_on_read_no_overlap(row, quality_threshold=quality_threshold, mode=mode),
            axis=1, result_type="expand")
    called_bases = pd.concat(list(bases_object[0]))
    ignored_bases = pd.concat(list(bases_object[1]))
    called_bases, ignored_bases = filter_multimapped_bases(called_bases, ignored_bases, reads_overlap)
    ignored_bases.to_csv(os.path.join(output_dir, base_filename+".ignored_bases"), sep="\t", index_label='ref_pos')
    called_bases.to_csv(os.path.join(output_dir, base_filename+".called_bases"), sep="\t", index_label='ref_pos')


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
    # TODO: test all these blast defaults.
    # Set blast defaults:
    if not task:
        task = "blastn"
    if evalue is None:
        evalue = 1e-7
    if not dust:
        dust = "no"
    if num_alignments is None:
        num_alignments = 1000000
    if not soft_masking:
        soft_masking = "F"
    if perc_identity is None:
        perc_identity = 85
    if not mode:
        mode = "RefToSeq"  # TODO: test differences between modes.
    run_blast(reads_fasta=reads_fasta_file_path, reference=reference, output=blast_output_file, mode=mode, task=task,
              evalue=evalue, num_alignments=num_alignments, dust=dust, soft_masking=soft_masking, log=log,
              perc_identity=perc_identity)
    basecall_output = os.path.join(output_dir, 'basecall')
    if not quality_threshold:
        quality_threshold = 30
    if reads_overlap == 'N':
        reads_overlap = False
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
    parser.add_argument("-q", "--quality_threshold", type=int,
                        help="bases with quality under this will be ignored (default: 30)")
    parser.add_argument("-bt", "--blast_task", help="blast's task parameter (default: blastn")
    parser.add_argument("-be", "--blast_evalue", help="blast's evalue parameter (default: 1e-7)", type=float)
    parser.add_argument("-bd", "--blast_dust", help="blast's dust parameter (default: on)")
    parser.add_argument("-bn", "--blast_num_alignments", type=int,
                        help="blast's num_alignments parameter (default: 1000000)")
    parser.add_argument("-bp", "--blast_perc_identity", type=int, help="blast's perc_identity parameter (default: 85)")
    parser.add_argument("-bs", "--blast_soft_masking", type=int, help="blast's soft_masking parameter (default: F)")
    parser.add_argument("-bm", "--blast_mode", help="RefToSeq or SeqToRef (default: RefToSeq)")  # TODO: docs
    args = parser.parse_args()
    process_fastq(fastq_file=args.fastq_file, reference=args.reference_path, output_dir=args.output_dir,
                  reads_overlap=args.reads_overlap, quality_threshold=args.quality_threshold, task=args.blast_task,
                  evalue=args.blast_evalue, dust=args.blast_dust, num_alignments=args.blast_num_alignments,
                  mode=args.blast_mode, perc_identity=args.blast_perc_identity, soft_masking=args.blast_soft_masking)
