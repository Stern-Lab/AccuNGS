import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def _get_total_read_list(mutations_reads_list):
    complete_read_list = set()
    for read_list in mutations_reads_list.read_id:
        complete_read_list |= eval(read_list)
    return complete_read_list


def get_mutations_linked_with_position(x, variants_list, mutation_read_list):
    ret = []
    relevant_positions = set(range(x+1, x+300)) & set(variants_list.get_level_values(0))
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


def get_variants_list(freqs_file, min_base_count):
    freqs = pd.read_csv(freqs_file, sep='\t')
    variants = freqs[(freqs.read_base != freqs.ref_base) & (freqs.base_count > min_base_count)].set_index(
        ['ref_pos', 'read_base'])
    return variants.index


def calculate_linked_mutations(freqs_file_path, output, mutation_read_list_path):
    # TODO: mp support if needed.
    variants_list = get_variants_list(freqs_file_path, min_base_count=3)  # TODO: set defaults.
    mutation_read_list = pd.read_csv(mutation_read_list_path, sep="\t").set_index(
        ['ref_pos', 'read_base'], drop=True)
    linked_mutations = {}
    for position in variants_list.get_level_values(0).astype(int):
        linked_mutations[position] = get_mutations_linked_with_position(position, variants_list=variants_list,
                                                                        mutation_read_list=mutation_read_list)
    pd.concat(linked_mutations.values()).to_csv(output, sep='\t', index=False)
