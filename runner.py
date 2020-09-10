import argparse
import os
import multiprocessing as mp
import pandas as pd
from Bio import pairwise2

from data_preparation import prepare_data
from graph_haplotypes import graph_haplotypes
from mutations_linking import get_variants_list, get_mutations_linked_with_position
from plotting import graph_summary
from processing import process_fastq
from aggregation import aggregate_processed_output, create_freqs_file
from logger import pipeline_logger
from utils import get_files_in_dir, get_sequence_from_fasta, get_mp_results_and_report, create_new_ref_with_freqs, \
    get_files_by_extension, concatenate_files_by_extension
from haplotypes.co_occurs_to_stretches import calculate_stretches


def parallel_process(processing_dir, fastq_files, reference_file, quality_threshold, task, evalue, dust, num_alignments,
                     soft_masking, perc_identity, mode, cpu_count):
    with mp.Pool(cpu_count) as pool:
        parts = [pool.apply_async(process_fastq, args=(fastq_file, reference_file, processing_dir, quality_threshold, task,
                                                       evalue, dust, num_alignments, soft_masking, perc_identity, mode,))
                 for fastq_file in fastq_files]
        get_mp_results_and_report(parts)


def get_stages_list(stages_range):
    stages_dict = {1: 'prepare data', 2: 'process data', 3: 'aggregate output', 4: 'compute haplotypes',
                   5: 'graph haplotypes'}
    if len(stages_range) == 2:
        stages = range(stages_range[0], stages_range[1] + 1)
    elif len(stages_range) == 1:
        stages = stages_range
    else:
        raise Exception("Parameter stages_range must be 1 or 2 numbers.")
    return [stages_dict[stage] for stage in stages]


def set_filenames(output_dir):
    filenames = {"freqs_file_path": os.path.join(output_dir, 'freqs.tsv'),
                 "linked_mutations_path": os.path.join(output_dir, 'linked_mutations.tsv'),
                 "mutation_read_list_path": os.path.join(output_dir, 'mutation_read_list.tsv'),
                 "stretches": os.path.join(output_dir, 'stretches.tsv'),
                 "blast_file": os.path.join(output_dir, 'blast.tsv'),
                 "read_counter_file": os.path.join(output_dir, 'read_counter.tsv'),
                 "summary_graphs": os.path.join(output_dir, 'summary.png')}
    return filenames


def calculate_linked_mutations(freqs_file_path, mutation_read_list, max_read_length, output_dir):
    variants_list = get_variants_list(freqs_file_path)
    for position in mutation_read_list.index.get_level_values(0).astype(int).unique():
        output_path = os.path.join(output_dir, f"{position}_linked_mutations.tsv")
        get_mutations_linked_with_position(position, variants_list=variants_list,
                                           mutation_read_list=mutation_read_list,
                                           max_read_size=max_read_length,
                                           output_path=output_path)


def parallel_calc_linked_mutations(freqs_file_path, output_dir, mutation_read_list_path, max_read_length, part_size,
                                   cpu_count):
    mutation_read_list = pd.read_csv(mutation_read_list_path, sep="\t").set_index(['ref_pos', 'read_base'], drop=True)
    positions = mutation_read_list.index.get_level_values(0).astype(int).unique()
    if max_read_length is None:
        max_read_length = 350
    if not part_size:
        part_size = 50  # smaller parts take less time to compute but take more time to aggregate and use more RAM.
    mutation_read_list_parts = {}
    start_index = 0
    end_index = 0
    while end_index < len(positions):
        next_start_index = start_index + part_size
        next_end_index = min(next_start_index + part_size + max_read_length, len(positions))
        if next_end_index == len(positions):
            end_index = next_end_index
        else:
            end_index = start_index + part_size + max_read_length
        start_position = positions[start_index]
        end_position = positions[end_index-1]
        mutation_read_list_parts[f"{start_index}_{end_index}"] = mutation_read_list.loc[start_position:end_position]
        start_index += part_size
    with mp.Pool(cpu_count) as pool:
        parts = [pool.apply_async(calculate_linked_mutations,
                                  args=(freqs_file_path, read_list, max_read_length, output_dir))
                 for read_list in mutation_read_list_parts.values()]
        get_mp_results_and_report(parts)


def check_consensus_alignment_with_ref(reference_file, with_indels, min_coverage, iteration_data_dir, basecall_dir,
                                       iteration_counter):
    reference = get_sequence_from_fasta(reference_file)
    freqs_file_path = os.path.join(iteration_data_dir, f"freqs_{iteration_counter}.tsv")
    called_bases_files = get_files_by_extension(basecall_dir, "called_bases")
    create_freqs_file(called_bases_files=called_bases_files, output_path=freqs_file_path)
    if with_indels == "Y":
        consensus_path = os.path.join(iteration_data_dir, f"consensus_with_indels_{iteration_counter}.fasta")
        drop_indels = False
    else:
        consensus_path = os.path.join(iteration_data_dir, f"consensus_without_indels_{iteration_counter}.fasta")
        drop_indels = True
    create_new_ref_with_freqs(reference_fasta_file=reference_file, freqs_file=freqs_file_path,
                              min_coverage=min_coverage, output_file=consensus_path, drop_indels=drop_indels)
    consensus = get_sequence_from_fasta(consensus_path)
    alignment_score = pairwise2.align.globalxx(consensus, reference, score_only=True)
    alignment_score = alignment_score / max(len(consensus), len(reference))
    return alignment_score


def get_consensus_path(basecall_iteration_counter, consolidate_consensus_with_indels, iteration_data_dir):
    consensus = os.path.join(iteration_data_dir, "consensus_with")
    if consolidate_consensus_with_indels != "Y":
        consensus += "out"
    consensus += f"_indels_{basecall_iteration_counter}.fasta"
    return consensus


def runner(input_dir, reference_file, output_dir, stages_range, max_basecall_iterations, part_size, min_coverage,
           quality_threshold, task, evalue, dust, num_alignments, soft_masking, perc_identity, mode, max_read_size,
           consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot, cleanup,
           cpu_count):

    # TODO: trim read_ids to save ram
    if not cpu_count:
        cpu_count = mp.cpu_count()
    os.makedirs(output_dir, exist_ok=True)
    log = pipeline_logger(logger_name='AccuNGS-Runner', log_folder=output_dir)
    log.debug(f"runner params: {locals()}")
    if max_basecall_iterations is None:
        max_basecall_iterations = 2   # TODO: what should this be?!
    filenames = set_filenames(output_dir)
    stages = get_stages_list(stages_range)
    log.info(f"Running stages: {stages}")
    processing_dir = os.path.join(output_dir, "processing")
    linked_mutations_dir = os.path.join(output_dir, "linked_mutations")
    if 'prepare data' in stages:
        data_dir = os.path.join(output_dir, "data")
        os.makedirs(data_dir, exist_ok=True)
        log.info("Preparing data")
        prepare_data(input_dir=input_dir, output_dir=data_dir, part_size=part_size, opposing_strings=None)
    else:
        data_dir = input_dir
    if 'process data' in stages:
        os.makedirs(processing_dir, exist_ok=True)
        data_files = get_files_in_dir(data_dir)
        fastq_files = [file_path for file_path in data_files if "fastq.part_" in os.path.basename(file_path)]
        log.info(f"Processing {len(fastq_files)} fastq files.")
        for basecall_iteration_counter in range(1, max_basecall_iterations + 1):
            log.info(f"Processing fastq files iteration {basecall_iteration_counter}/{max_basecall_iterations}")
            # TODO: whats up with the different modes?!?
            parallel_process(processing_dir=processing_dir, fastq_files=fastq_files, reference_file=reference_file,
                             quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust, mode=mode,
                             num_alignments=num_alignments, soft_masking=soft_masking, perc_identity=perc_identity,
                             cpu_count=cpu_count)
            iteration_data_dir = os.path.join(output_dir, 'iteration_data')
            os.makedirs(iteration_data_dir, exist_ok=True)
            alignment_score = check_consensus_alignment_with_ref(reference_file=reference_file,
                                                                 iteration_counter=basecall_iteration_counter,
                                                                 basecall_dir=os.path.join(processing_dir, 'basecall'),
                                                                 with_indels=consolidate_consensus_with_indels,
                                                                 iteration_data_dir=iteration_data_dir,
                                                                 min_coverage=min_coverage)
            log.info(f'Iteration alignment score: {round(alignment_score,4)}')
            if alignment_score == 1:
                break
            consensus_path = get_consensus_path(basecall_iteration_counter, consolidate_consensus_with_indels,
                                                iteration_data_dir)
            reference_file = consensus_path
    if 'aggregate output' in stages:
        log.info("Aggregating processed fastq files outputs...")
        aggregate_processed_output(input_dir=processing_dir, output_dir=output_dir, cleanup=cleanup,
                                   reference=reference_file, min_coverage=min_coverage)
        graph_summary(freqs_file=filenames['freqs_file_path'], blast_file=filenames['blast_file'],
                      read_counter_file=filenames['read_counter_file'], stretches_file=filenames['stretches'],
                      output_file=filenames['summary_graphs'], min_coverage=min_coverage,
                      stretches_to_plot=stretches_to_plot)  # TODO: drop low quality mutations?
        log.info(f"Most outputs are ready in {output_dir} !")
    if 'compute haplotypes' in stages:
        log.info(f"Calculating linked mutations...")
        os.makedirs(linked_mutations_dir, exist_ok=True)
        # TODO: optimize part size
        parallel_calc_linked_mutations(freqs_file_path=filenames['freqs_file_path'], cpu_count=cpu_count,
                                       mutation_read_list_path=filenames['mutation_read_list_path'],
                                       output_dir=linked_mutations_dir, max_read_length=max_read_size,
                                       part_size=100)  # TODO: drop low quality mutations?, set part_size as param.
    if 'graph haplotypes' in stages:
        log.info(f"Aggregating linked mutations to stretches...")
        concatenate_files_by_extension(input_dir=linked_mutations_dir, extension='tsv',
                                       output_path=filenames['linked_mutations_path'])
        calculate_stretches(filenames['linked_mutations_path'], max_pval=stretches_pvalue, distance=stretches_distance,
                            output=filenames['stretches'])  #TODO: refactor that function
        graph_summary(freqs_file=filenames['freqs_file_path'], blast_file=filenames['blast_file'],
                      read_counter_file=filenames['read_counter_file'], stretches_file=filenames['stretches'],
                      output_file=filenames['summary_graphs'], min_coverage=min_coverage,
                      stretches_to_plot=stretches_to_plot)  # TODO: drop low quality mutations?
        graph_haplotypes(input_file=filenames['stretches'], number_of_stretches=stretches_to_plot,
                         output_dir=output_dir)
        log.info(f"Done! Output files are in directory: {output_dir}")

    #TODO: test everything, finish up, pbs_runner?


def create_runner_parser():
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
    parser.add_argument("-mc", "--min_coverage", type=int, default=10,
                        help="minimal coverage to plot in summary graphs (default: 10)")  # TODO: different default?
    parser.add_argument("-ccwi", "--consolidate_consensus_with_indels", type=str, default="Y",
                        help="Y/N where N means we consolidate consensus without indels (default: Y)")
    parser.add_argument("-sp", "--stretches_pvalue", type=float,
                        help="only consider joint mutations with pvalue below this (default: 10**-9")  # TODO: better docs
    parser.add_argument("-sd", "--stretches_distance", type=float,
                        help="mean transitive distance between joint mutations to calculate stretches (default: 10)")
    parser.add_argument("-stp", "--stretches_to_plot", type=int, default=5,
                        help="number of stretches to plot in deep dive (default: 5)")
    parser.add_argument("-smrs", "--stretches_max_read_size", type=int,
                        help="look this many positions forward for joint mutations (default: 350)")
    parser.add_argument("-c", "--cleanup", help="remove input folder when done (default: Y)", default="Y")
    parser.add_argument("-cc", "--cpu_count", help="max number of cpus to use (default: all)", type=int)
    return parser


if __name__ == "__main__":
    parser = create_runner_parser()
    args = parser.parse_args()
    runner(input_dir=args.input_dir, output_dir=args.output_dir, reference_file=args.reference_file,
           stages_range=args.stages_range, max_basecall_iterations=args.max_basecall_iterations,
           part_size=args.part_size, quality_threshold=args.quality_threshold, task=args.blast_task,
           evalue=args.blast_evalue, dust=args.blast_dust, num_alignments=args.blast_num_alignments,
           mode=args.blast_mode, perc_identity=args.blast_perc_identity, soft_masking=args.blast_soft_masking,
           min_coverage=args.min_coverage, consolidate_consensus_with_indels=args.consolidate_consensus_with_indels,
           stretches_pvalue=args.stretches_pvalue, stretches_distance=args.stretches_distance, cleanup=args.cleanup,
           stretches_to_plot=args.stretches_to_plot, max_read_size=args.stretches_max_read_size,
           cpu_count=args.cpu_count)
