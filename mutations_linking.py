import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from utils import drange

# TODO: combine with co-occurs_to_stretches.py
#       docs


def _get_total_read_list(mutations_reads_list):
    complete_read_list = set()
    for read_list in mutations_reads_list:
        complete_read_list |= read_list
    return complete_read_list


def get_mutations_linked_with_position(x, variants_list, mutation_read_list, max_read_size, output_path):
    # TODO: speed!
    # TODO: test insertions
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    if max_read_size is None:
        max_read_size = 350
    ret = []
    x = round(x)
    relevant_positions = set(drange(x+1, x+max_read_size+0.1, 0.1)) & \
                         set(variants_list.get_level_values(0).astype(float)) & \
                         set(mutation_read_list.index.get_level_values(0).astype(float))
    mutation_read_list = mutation_read_list.read_id.map(eval)  # turn all those strings into sets TODO: security?
    x_and_insertions = [pos for pos in variants_list.get_level_values(0).astype(float) if (pos < x+1) and (pos >= x)]
    for x in x_and_insertions:
        x_read_list = _get_total_read_list(mutation_read_list.loc[x])
        for x_mutation in mutation_read_list.loc[x].index:
            if (x, x_mutation) in variants_list:
                x_mutation_read_list = mutation_read_list.loc[x, x_mutation]
                for y in relevant_positions:
                    y_read_list = _get_total_read_list(mutation_read_list.loc[y])
                    for y_mutation in mutation_read_list.loc[y].index:
                        y_mutation_read_list = mutation_read_list.loc[(y, y_mutation)]
                        both_mutations_on_read = len(set(x_mutation_read_list) & set(y_mutation_read_list))
                        if both_mutations_on_read > 0:
                            just_x_on_read = len(set(x_mutation_read_list) - set(y_mutation_read_list))
                            just_y_on_read = len(set(y_mutation_read_list) - set(x_mutation_read_list))
                            both_mutations_not_on_read = len((x_read_list | y_read_list) - (
                                set(x_mutation_read_list) | set(y_mutation_read_list)))
                            contingency_table = [[both_mutations_on_read, just_x_on_read],
                                                 [just_y_on_read, both_mutations_not_on_read]]
                            p_value = fisher_exact(contingency_table, alternative='greater')[1]
                            frequency = both_mutations_on_read / np.sum(contingency_table)
                            ret.append({'x_pos': x,
                                        'x_mutation': x_mutation,
                                        'y_pos': y,
                                        'y_mutation': y_mutation,
                                        'p_value': p_value,
                                        'frequency': frequency})
    pd.DataFrame.from_dict(ret).to_csv(output_path, sep='\t', index=False)


def get_variants_list(freqs_file):
    freqs = pd.read_csv(freqs_file, sep='\t')
    freqs['ref_pos'] = freqs['ref_pos'].round(3)
    variants = freqs[(freqs['base_rank'] != 0) & (freqs.base_count > 0)].set_index(['ref_pos', 'read_base'])
    return variants.index


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--position", type=float, help="position to check joint correlations with")
    parser.add_argument("-f", "--freqs_file", help="freqs file path to get variant list from")
    parser.add_argument("-r", "--mutation_read_list", help="mutation read list file path")
    parser.add_argument("-m", "--max_read_size", help="length of largest read - will look this many bases forward",
                        type=int)
    parser.add_argument("-o", "--output_path", help="where the output file goes")
    args = parser.parse_args()
    variants_list = get_variants_list(args.freqs_file)
    mutation_read_list = pd.read_csv(args.mutation_read_list, sep="\t").set_index(['ref_pos', 'read_base'], drop=True)
    get_mutations_linked_with_position(x=args.position, max_read_size=args.max_read_size, output_path=args.output_path,
                                       mutation_read_list=mutation_read_list, variants_list=variants_list)
