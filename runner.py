import argparse
import os
import multiprocessing as mp
import pandas as pd
from Bio import pairwise2

from data_preparation import prepare_data
from graph_haplotypes import graph_haplotypes
from haplotypes.mutations_linking import get_variants_list, get_mutations_linked_with_position
from plotting import graph_summary
from processing import process_fastq
from aggregation import aggregate_processed_output
from logger import pipeline_logger
from utils import get_files_in_dir, get_sequence_from_fasta, get_mp_results_and_report
from haplotypes.co_occurs_to_stretches import calculate_stretches


def parallel_process(processing_dir, fastq_files, reference_file, quality_threshold, task, evalue, dust, num_alignments,
                     soft_masking, perc_identity, mode):
    pool = mp.Pool(processes=mp.cpu_count())
    parts = [pool.apply_async(process_fastq, args=(fastq_file, reference_file, processing_dir, quality_threshold, task,
                                                   evalue, dust, num_alignments, soft_masking, perc_identity, mode,))
             for fastq_file in fastq_files]
    get_mp_results_and_report(parts)


def get_stages_list(stages_range):
    stages_dict = {1: 'prepare_data', 2: 'process_data', 3: 'analyze'}
    if len(stages_range) == 2:
        stages = range(stages_range[0], stages_range[1] + 1)
    elif len(stages_range) == 1:
        stages = stages_range
    else:
        raise Exception("Parameter stages_range must be 1 or 2 numbers.")
    return [stages_dict[stage] for stage in stages]


def get_alignment_score(consensus, reference_file):
    reference = get_sequence_from_fasta(reference_file)
    alignment_score = pairwise2.align.globalxx(consensus, reference, score_only=True)
    alignment_score = alignment_score / max(len(consensus), len(reference))
    return alignment_score


def set_filenames(output_dir):
    filenames = {"freqs_file_path": os.path.join(output_dir, 'freqs.tsv'),
                 "linked_mutations_path": os.path.join(output_dir, 'linked_mutations.tsv'),
                 "mutation_read_list_path": os.path.join(output_dir, 'mutation_read_list.tsv'),
                 "stretches": os.path.join(output_dir, 'streches.tsv'),
                 "blast_file": os.path.join(output_dir, 'blast.tsv'),
                 "read_counter_file": os.path.join(output_dir, 'read_counter.tsv'),
                 "summary_graphs": os.path.join(output_dir, 'summary.png')}
    return filenames


def calculate_linked_mutations(freqs_file_path, mutation_read_list, max_read_length):
    variants_list = get_variants_list(freqs_file_path)
    linked_mutations = []
    for position in mutation_read_list.index.get_level_values(0).astype(int).unique():
        linked_mutations.append(get_mutations_linked_with_position(position, variants_list=variants_list,
                                                                   mutation_read_list=mutation_read_list,
                                                                   max_read_size=max_read_length))
    return pd.concat(linked_mutations)


def parallel_calc_linked_mutations(freqs_file_path, output, mutation_read_list_path, max_read_length, part_size):
    mutation_read_list = pd.read_csv(mutation_read_list_path, sep="\t").set_index(['ref_pos', 'read_base'], drop=True)
    print(mutation_read_list)
    positions = mutation_read_list.index.get_level_values(0).astype(int).unique()
    if max_read_length is None:
        max_read_length = 350
    if not part_size:
        part_size = 50  # smaller parts take less time to compute but take more time to aggregate and use more RAM.
    mutation_read_list_parts = {}
    start_index = min(positions)
    end_index = 0

    while end_index < max(positions):
        next_start_index = start_index + part_size
        next_end_index = min(next_start_index + part_size + max_read_length, max(positions))
        if next_end_index == max(positions):
            end_index = next_end_index
        else:
            end_index = start_index + part_size + max_read_length
        print(start_index, end_index)
        mutation_read_list_parts[f"{start_index}_{end_index}"] = mutation_read_list.loc[start_index:end_index]
        start_index += part_size
    pool = mp.Pool(processes=mp.cpu_count())
    parts = [pool.apply_async(calculate_linked_mutations,
                              args=(freqs_file_path, read_list, max_read_length))
             for read_list in mutation_read_list_parts.values()]
    results = get_mp_results_and_report(parts)
    pd.concat(results).to_csv(output, sep='\t', index=False)


# TODO: trim read_ids to save ram

def runner(input_dir, reference_file, output_dir, stages_range, max_basecall_iterations, part_size, min_coverage,
           quality_threshold, task, evalue, dust, num_alignments, soft_masking, perc_identity, mode, max_read_size,
           consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot):
    log = pipeline_logger(logger_name='AccuNGS-Runner', log_folder=output_dir)
    log.debug(f"runner params: {locals()}")
    if max_basecall_iterations is None:
        max_basecall_iterations = 2   # TODO: what should this be?!
    filenames = set_filenames(output_dir)
    stages = get_stages_list(stages_range)
    log.info(f"Running stages: {stages}")
    if 'prepare_data' in stages:
        data_dir = os.path.join(output_dir, "data")
        os.makedirs(data_dir, exist_ok=True)
        log.info("Preparing data")
        prepare_data(input_dir=input_dir, output_dir=data_dir, part_size=part_size, opposing_strings=None)
    else:
        data_dir = input_dir
    if 'process_data' in stages:
        processing_dir = os.path.join(output_dir, "processing")
        os.makedirs(processing_dir, exist_ok=True)
        data_files = get_files_in_dir(data_dir)
        fastq_files = [file_path for file_path in data_files if "fastq.part_" in os.path.basename(file_path)]
        if consolidate_consensus_with_indels == "Y":
            consensus_file = os.path.join(output_dir, 'consensus_with_indels.fasta')
        else:
            consensus_file = os.path.join(output_dir, 'consensus_without_indels.fasta')
        log.info(f"Processing {len(fastq_files)} fastq files.")
        for basecall_iteration_counter in range(1, max_basecall_iterations + 1):
            #TODO: test that it sometimes takes more than just 2 iterations.
            log.info(f"Processing fastq files iteration {basecall_iteration_counter}/{max_basecall_iterations}")
            # TODO: whats up with the different modes?!?
            parallel_process(processing_dir=processing_dir, fastq_files=fastq_files, reference_file=reference_file,
                             quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust, mode=mode,
                             num_alignments=num_alignments, soft_masking=soft_masking, perc_identity=perc_identity)
            log.info("Aggregating processed fastq files outputs...")
            aggregate_processed_output(input_dir=processing_dir, output_dir=output_dir,
                                       reference=reference_file, min_coverage=min_coverage)
            # TODO: factor out creating consensus from aggregation
            consensus = get_sequence_from_fasta(consensus_file)
            alignment_score = get_alignment_score(consensus=consensus, reference_file=reference_file)
            log.info(f"Iteration done with alignment score: {round(alignment_score,4)}")
            if alignment_score == 1:
                break
            reference_file = consensus_file
        log.info(f"Calculating linked mutations...")
        parallel_calc_linked_mutations(freqs_file_path=filenames['freqs_file_path'],
                                       mutation_read_list_path=filenames['mutation_read_list_path'],
                                       output=filenames['linked_mutations_path'], max_read_length=max_read_size,
                                       part_size=30)  # TODO: drop low quality mutations?, set part_size as param.
        log.info(f"Aggregating linked mutations to stretches...")
        calculate_stretches(filenames['linked_mutations_path'], max_pval=stretches_pvalue, distance=stretches_distance,
                            output=filenames['stretches'])  #TODO: refactor that function
        log.info(f"Processing done. Output files are in directory: {output_dir}")
    if 'analyze' in stages:
        log.info("Drawing summary graphs!")
        graph_summary(freqs_file=filenames['freqs_file_path'], blast_file=filenames['blast_file'],
                      read_counter_file=filenames['read_counter_file'], stretches_file=filenames['stretches'],
                      output_file=filenames['summary_graphs'], min_coverage=min_coverage,
                      stretches_to_plot=stretches_to_plot)  # TODO: drop low quality mutations?
        graph_haplotypes(input_file=filenames['stretches'], number_of_stretches=stretches_to_plot,
                         output_dir=output_dir)
    #TODO: test everything, finish up, pbs_runner?
    log.info("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing basecall files")
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-r", "--reference_file", required=True)
    parser.add_argument("-s", "--stages_range", nargs="+", type=int, help="start and end stages separated by spaces")
    parser.add_argument("-m", "--max_basecall_iterations", type=int,
                        help="number of times to run basecall before giving up equalizing reference with consensus")
    parser.add_argument("-p", "--part_size", type=int,
                        help="size of part of data file in string length")
    parser.add_argument("-bt", "--blast_task", help="blast's task parameter (default: blastn")
    parser.add_argument("-be", "--blast_evalue", help="blast's evalue parameter (default: 1e-7)", type=float)
    parser.add_argument("-bd", "--blast_dust", help="blast's dust parameter (default: on)")
    parser.add_argument("-bn", "--blast_num_alignments", type=int,
                        help="blast's num_alignments parameter (default: 1000000)")
    parser.add_argument("-bp", "--blast_perc_identity", type=int, help="blast's perc_identity parameter (default: 85)")
    parser.add_argument("-bs", "--blast_soft_masking", type=int, help="blast's soft_masking parameter (default: F)")
    parser.add_argument("-bm", "--blast_mode", help="RefToSeq or SeqToRef (default: RefToSeq)")  # TODO: docs
    parser.add_argument("-qt", "--quality_threshold", type=int,
                        help="phred score must be higher than this to be included (default: 30)")
    parser.add_argument("-mc", "--min_coverage", type=int,
                        help="minimal coverage to plot in summary graphs (default: 10)") # TODO: different default?
    parser.add_argument("-ccwi", "--consolidate_consensus_with_indels", type=str, default="Y",
                        help="Y/N where N means we consolidate consensus without indels (default: Y)")
    parser.add_argument("-sp", "--stretches_pvalue", type=float,
                        help="only consider joint mutations with pvalue below this (default: 10**-9")  #TODO: better docs
    parser.add_argument("-sd", "--stretches_distance", type=float,
                        help="mean transitive distance between joint mutations to calculate stretches (default: 10)")
    parser.add_argument("-stp", "--stretches_to_plot", type=int, default=5,
                        help="number of stretches to plot in deep dive (default: 5)")
    parser.add_argument("-smrs", "--stretches_max_read_size", type=int,
                        help="look this many positions forward for joint mutations (default: 350)")

    args = parser.parse_args()
    runner(input_dir=args.input_dir, output_dir=args.output_dir, reference_file=args.reference_file,
           stages_range=args.stages_range, max_basecall_iterations=args.max_basecall_iterations,
           part_size=args.part_size, quality_threshold=args.quality_threshold, task=args.blast_task,
           evalue=args.blast_evalue, dust=args.blast_dust, num_alignments=args.blast_num_alignments,
           mode=args.blast_mode, perc_identity=args.blast_perc_identity, soft_masking=args.blast_soft_masking,
           min_coverage=args.min_coverage, consolidate_consensus_with_indels=args.consolidate_consensus_with_indels,
           stretches_pvalue=args.stretches_pvalue, stretches_distance=args.stretches_distance,
           stretches_to_plot=args.stretches_to_plot, max_read_size=args.stretches_max_read_size)
