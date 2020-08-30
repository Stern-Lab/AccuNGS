import argparse
import os
import numpy as np
import pandas as pd

from utils import get_files_by_extension, concatenate_files_by_extension, create_new_ref_with_freqs


def convert_called_bases_to_freqs(called_bases):
    dummy_bases = []
    for pos in called_bases.ref_pos.unique():
        for base in ['A', 'G', 'T', 'C', '-']:
            dummy_bases.append({'ref_pos': pos, 'read_base': base})
    freq_dummies = pd.DataFrame.from_dict(dummy_bases)
    freqs = pd.concat([freq_dummies,called_bases])
    freqs = freqs.groupby('ref_pos').read_base.value_counts()-1
    ref_df = called_bases[['ref_pos', 'ref_base']]
    return freqs, ref_df


def aggregate_called_bases(called_bases_files):
    freqs = pd.Series(dtype=int)
    ref_df = pd.DataFrame()
    for called_bases_file in called_bases_files:
        called_bases_df = pd.read_csv(called_bases_file, sep="\t")
        freqs_part, ref_df_part = convert_called_bases_to_freqs(called_bases_df)
        if freqs.empty:
            freqs = freqs_part
        else:
            freqs = freqs.add(freqs_part, fill_value=0)
            ref_df = pd.concat([ref_df, ref_df_part]).drop_duplicates()
    freqs.name = 'base_count'
    freqs = pd.DataFrame(freqs).reset_index()
    freqs = freqs.merge(ref_df, on=['ref_pos'], how='left')
    return freqs


def create_freqs_file(called_bases_files, output_path):
    freqs = aggregate_called_bases(called_bases_files)
    coverage = freqs.groupby('ref_pos').base_count.sum()
    freqs['coverage'] = freqs.ref_pos.map(lambda pos: coverage[pos])
    freqs['frequency'] = freqs['base_count'] / freqs['coverage']
    freqs['base_rank'] = 5 - freqs.groupby('ref_pos').base_count.rank('min')
    freqs.to_csv(output_path, sep="\t", index=False)


def collect_reads_from_row(row):
    read_lists = {'other_reads': row.read_id, 'these_reads': row.read_id_this}
    for key, read_list in read_lists.items():
        if read_list is np.nan:
            read_lists[key] = ""
        else:
            read_lists[key] = read_list
    return set(read_lists['other_reads']) | set(read_lists['these_reads'])


def append_read_lists(read_list, this_read_list):
    read_list = read_list.join(this_read_list, rsuffix="_this", how='outer')
    read_list['read_id'] = read_list.apply(collect_reads_from_row, axis=1)
    return read_list[['read_id']]


def create_mutation_read_list_file(called_bases_files, output_path):
    read_list = pd.DataFrame()
    for bases_file in called_bases_files:
        this_read_list = pd.read_csv(bases_file, sep="\t").groupby(['ref_pos', 'read_base']).read_id.unique()
        if read_list.empty:
            read_list = pd.DataFrame(this_read_list)
        else:
            read_list = append_read_lists(read_list=read_list, this_read_list=this_read_list)
    read_list.to_csv(output_path, sep='\t')


def aggregate_computation_output(input_dir, output_dir, reference, min_coverage=None):
    # TODO: graphs.
    if not min_coverage:
        min_coverage = 10
    os.makedirs(output_dir, exist_ok=True)
    freqs_file_path = os.path.join(output_dir, "freqs.tsv")
    mutation_read_list_path = os.path.join(output_dir, "mutation_read_list.tsv")
    basecall_dir = os.path.join(input_dir, 'basecall')
    blast_dir = os.path.join(input_dir, 'blast')
    called_bases_files = get_files_by_extension(basecall_dir, "called_bases")
    if len(called_bases_files) == 0:
        raise Exception(f"Could not find files of type *.called_bases in {input_dir}")
    create_freqs_file(called_bases_files=called_bases_files, output_path=freqs_file_path)
    create_mutation_read_list_file(called_bases_files=called_bases_files, output_path=mutation_read_list_path)
    concatenate_files_by_extension(input_dir=blast_dir, extension="blast",
                                   output_path=os.path.join(output_dir, "blast.tsv"))
    concatenate_files_by_extension(input_dir=basecall_dir, extension="called_bases",
                                   output_path=os.path.join(output_dir, "called_bases.tsv"))
    concatenate_files_by_extension(input_dir=basecall_dir, extension="ignored_bases",
                                   output_path=os.path.join(output_dir, "ignored_bases.tsv"))
    concatenate_files_by_extension(input_dir=basecall_dir, extension="suspicious_reads",
                                   output_path=os.path.join(output_dir, "suspicious_reads.tsv"))
    concatenate_files_by_extension(input_dir=basecall_dir, extension="ignored_reads",
                                   output_path=os.path.join(output_dir, "ignored_reads.tsv"))
    concatenate_files_by_extension(input_dir=basecall_dir, extension="ignored_reads",
                                   output_path=os.path.join(output_dir, "ignored_reads.tsv"))
    create_new_ref_with_freqs(reference_fasta_file=reference, freqs_file=freqs_file_path, min_coverage=min_coverage,
                              output_file=os.path.join(output_dir, "consensus_without_indels.fasta"), drop_indels=True)
    create_new_ref_with_freqs(reference_fasta_file=reference, freqs_file=freqs_file_path, min_coverage=min_coverage,
                              output_file=os.path.join(output_dir, "consensus_with_indels.fasta"), drop_indels=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing basecall files")
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-r", "--reference_file", required=True)
    parser.add_argument("-mc", "--min_coverage",
                        help="bases with less than this coverage will be excluded from affecting the consensus "
                             "(default: 10)")
    args = parser.parse_args()
    aggregate_computation_output(input_dir=args.input_dir, output_dir=args.output_dir, reference=args.reference,
                                 min_coverage=args.min_coverage)
