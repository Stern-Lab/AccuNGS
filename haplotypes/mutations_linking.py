import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def _get_total_read_list(mutations_reads_list):
    complete_read_list = set()
    for read_list in mutations_reads_list.read_id:
        complete_read_list |= eval(read_list)
    return complete_read_list


def get_mutations_linked_with_position(x, variants_list, mutation_read_list, max_read_size):
    if max_read_size is None:
        max_read_size = 350
    ret = []
    relevant_positions = set(range(x+1, x+max_read_size)) & set(variants_list.get_level_values(0).astype(int)) & set(
        mutation_read_list.index.get_level_values(0).astype(int))
    x_read_list = _get_total_read_list(mutation_read_list.loc[x])
    for x_mutation in mutation_read_list.loc[x].index:
        if (x, x_mutation) in variants_list:
            x_mutation_read_list = eval(mutation_read_list.loc[x, x_mutation].values[0])
            for y in relevant_positions:
                y_read_list = _get_total_read_list(mutation_read_list.loc[y])
                for y_mutation in mutation_read_list.loc[y].index:
                    y_mutation_read_list = eval(mutation_read_list.loc[(y, y_mutation)].values[0])
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
    return pd.DataFrame.from_dict(ret)


def get_variants_list(freqs_file):
    freqs = pd.read_csv(freqs_file, sep='\t')
    variants = freqs[(freqs['base_rank'] != 0) & (freqs.base_count > 0)].set_index(['ref_pos', 'read_base'])
    return variants.index
