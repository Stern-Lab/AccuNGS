"""
This script is used by runner.py in order to aggregate all the parallel computation output done by instances of
processing.py.

It trims read_ids (to save space), concats the blast and basecalling outputs, creates the mutation_read_list (called
separately by the runner) and creates the famous freqs file and the new reference files from it.
"""
import argparse
import json
import os
import numpy as np
import pandas as pd

from utils import get_files_by_extension, concatenate_files_by_extension, create_new_ref_with_freqs, get_files_in_dir


def convert_called_bases_to_freqs(called_bases):
    dummy_bases = []
    for pos in called_bases.ref_pos.unique():
        for base in ['A', 'G', 'T', 'C', '-']:
            dummy_bases.append({'ref_pos': pos, 'read_base': base})
    freq_dummies = pd.DataFrame.from_dict(dummy_bases)
    freqs = pd.concat([freq_dummies, called_bases])
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
    freqs['ref_pos'] = round(freqs['ref_pos'], 3)                         # fix that floating point nonsense
    return freqs


def create_freqs_file(called_bases_files, output_path):
    freqs = aggregate_called_bases(called_bases_files)
    freqs['round_ref_pos'] = round(freqs['ref_pos'])                   # to fix insertion coverage
    coverage = freqs.groupby('round_ref_pos').base_count.sum()
    freqs.drop(columns=['round_ref_pos'], inplace=True)
    freqs['coverage'] = freqs.ref_pos.map(lambda pos: coverage[round(pos)])
    freqs['frequency'] = freqs['base_count'] / freqs['coverage']
    freqs['base_rank'] = 5 - freqs.groupby('ref_pos').base_count.rank('min')
    freqs['probability'] = 1 - 10**(np.log10(1 - freqs["frequency"] + 1e-07) * (freqs["coverage"] + 1))
    # TODO: does probability logic make sense?
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
    read_list = read_list.reset_index()
    read_list['ref_pos'] = round(read_list['ref_pos'], 3)
    read_list.set_index(['ref_pos', 'read_base'], inplace=True)
    read_list.to_csv(output_path, sep='\t')


def aggregate_read_counters(read_counters, output_path):
    counter = {}
    for read_counter in read_counters:
        counter[read_counter] = pd.read_csv(read_counter, sep='\t')
    counters = pd.concat(counter.values())
    counters.groupby('read_id')['number_of_alignments'].sum().to_csv(output_path, sep='\t')


def update_prefix_dict(json_file, prefixes):
    if os.path.isfile(json_file):
        with open(json_file) as read_handle:
            read_id_prefix_dict = json.load(read_handle)
        next_prefix_value = 1
        if len(read_id_prefix_dict) > 0:
            next_prefix_value = max(read_id_prefix_dict.values()) + 1
    else:
        read_id_prefix_dict = {}
        next_prefix_value = 1
    for prefix in prefixes:
        if prefix not in read_id_prefix_dict.keys():
            read_id_prefix_dict[prefix] = next_prefix_value
            next_prefix_value += 1
    with open(json_file, 'w') as write_handle:
        json.dump(read_id_prefix_dict, write_handle)
    return read_id_prefix_dict


def trim_read_id_prefixes(files, read_id_prefix_file):
    prefix_length = 31
    for file in files:
        df = pd.read_table(file)
        if not df.empty:
            prefixes = df.read_id.str[:prefix_length].unique()
            prefix_dict = update_prefix_dict(read_id_prefix_file, prefixes)
    for file in files:
        df = pd.read_table(file)
        if not df.empty:
            df['read_id'] = df.read_id.map(lambda x: str(prefix_dict[x[:prefix_length]]) + x[prefix_length:])
            if file.endswith("bases"):
                df.to_csv(file, sep='\t', index=False, index_label='ref_pos')
            else:
                df.to_csv(file, sep='\t', index=False)


def aggregate_processed_output(input_dir, output_dir, reference, min_coverage):
    if not min_coverage:
        min_coverage = 10
    os.makedirs(output_dir, exist_ok=True)
    freqs_file_path = os.path.join(output_dir, "freqs.tsv")
    basecall_dir = os.path.join(input_dir, 'basecall')
    blast_dir = os.path.join(input_dir, 'blast')
    called_bases_files = get_files_by_extension(basecall_dir, "called_bases")
    if len(called_bases_files) == 0:
        raise Exception(f"Could not find files of type *.called_bases in {input_dir}")
    basecall_files = get_files_in_dir(basecall_dir)
    read_id_prefix_file = os.path.join(output_dir, "read_id_prefixes.json")
    trim_read_id_prefixes(files=basecall_files, read_id_prefix_file=read_id_prefix_file)
    create_freqs_file(called_bases_files=called_bases_files, output_path=freqs_file_path)
    read_counters = get_files_by_extension(basecall_dir, "read_counter")
    aggregate_read_counters(read_counters=read_counters, output_path=os.path.join(output_dir, "read_counter.tsv"))
    concatenate_files_by_extension(input_dir=blast_dir, extension="blast", remove_headers=False,
                                   output_path=os.path.join(output_dir, "blast.tsv"))
    for file_type in ['called_bases', 'ignored_bases', 'suspicious_reads', 'ignored_reads']:
        concatenate_files_by_extension(input_dir=basecall_dir, extension=file_type,
                                       output_path=os.path.join(output_dir, f"{file_type}.tsv"))
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
    aggregate_processed_output(input_dir=args.input_dir, output_dir=args.output_dir, reference=args.reference,
                               min_coverage=args.min_coverage)
