import argparse
import os
import pandas as pd

from utils import get_files_by_extension, concatenate_files_by_extension, make_reference_from_freqs


def convert_called_bases_to_freqs(called_bases):
    dummy_bases = []
    for pos in called_bases.ref_pos.unique():
        for base in ['A', 'G', 'T', 'C', '-']:
            dummy_bases.append({'ref_pos': pos, 'read_seq_plus': base})
    freq_dummies = pd.DataFrame.from_dict(dummy_bases)
    freqs = pd.concat([freq_dummies,called_bases]).rename(columns={'read_seq_plus': 'read_base'})
    freqs = freqs.groupby('ref_pos').read_base.value_counts()-1
    return freqs


def aggregate_called_bases(called_bases_files):
    freqs = pd.Series(dtype=int)
    for called_bases_file in called_bases_files:
        called_bases_df = pd.read_csv(called_bases_file, sep="\t")
        freqs_part = convert_called_bases_to_freqs(called_bases_df)
        if freqs.empty:
            freqs = freqs_part
        else:
            freqs = freqs.add(freqs_part, fill_value=0)
    return freqs


def create_freqs_file(input_dir, output_path):
    called_bases_files = get_files_by_extension(input_dir, "called_bases")
    if len(called_bases_files)==0:
        raise Exception(f"Could not find files of type *.called_bases in {input_dir}")
    freqs = aggregate_called_bases(called_bases_files)
    freqs.name = 'base_count'
    freqs = pd.DataFrame(freqs).reset_index()
    coverage = freqs.groupby('ref_pos').base_count.sum()
    freqs['coverage'] = freqs.ref_pos.map(lambda pos: coverage[pos])
    freqs['frequency'] = freqs['base_count'] / freqs['coverage']
    freqs['base_rank'] = 5 - freqs.groupby('ref_pos').base_count.rank('min')
    freqs.to_csv(output_path, sep="\t", index=False)


def aggregate_computation_output(input_dir, output_dir, reference, min_coverage=None):
    if not min_coverage:
        min_coverage = 10
    os.makedirs(output_dir, exist_ok=True)
    freqs_file_path = os.path.join(output_dir, "freqs.tsv")
    create_freqs_file(input_dir=input_dir, output_path=freqs_file_path)
    concatenate_files_by_extension(input_dir=input_dir, extension="called_bases",
                                   output_path=os.path.join(output_dir, "called_bases.tsv"))
    concatenate_files_by_extension(input_dir=input_dir, extension="ignored_bases",
                                   output_path=os.path.join(output_dir, "ignored_bases.tsv"))
    concatenate_files_by_extension(input_dir=input_dir, extension="suspicious_reads",
                                   output_path=os.path.join(output_dir, "suspicious_reads.tsv"))
    concatenate_files_by_extension(input_dir=input_dir, extension="ignored_reads",
                                   output_path=os.path.join(output_dir, "ignored_reads.tsv"))
    make_reference_from_freqs(reference_fasta_file=reference, freqs_file=freqs_file_path, min_coverage=min_coverage,
                              output_file=os.path.join(output_dir, "consensus_without_indels.fasta"), drop_indels=True)
    make_reference_from_freqs(reference_fasta_file=reference, freqs_file=freqs_file_path, min_coverage=min_coverage,
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
