import argparse
import os
import multiprocessing as mp
import sys

from Bio import pairwise2

from data_preparation import prepare_data
from graph_haplotypes import graph_haplotypes
from haplotypes.mutations_linking import calculate_linked_mutations
from plotting import graph_summary
from processing import compute
from aggregation import aggregate_computation_output
from logger import pipeline_logger
from utils import get_files_in_dir, get_sequence_from_fasta
from haplotypes.co_occurs_to_stretches import calculate_stretches


def parallel_compute(computation_dir, fastq_files, reference, max_cpus=None):
    if not max_cpus:
        max_cpus = mp.cpu_count()
    pool = mp.Pool(processes=max_cpus)
    parts = [pool.apply_async(compute, args=(fastq_file, reference, computation_dir, 30, False))
             for fastq_file in fastq_files]
    i = 1
    for res in parts:
        res.get()
        sys.stdout.write(f"\rDone {i}/{len(parts)} parts.")
        sys.stdout.flush()
        i = i + 1
    sys.stdout.write("\n")


def get_stages_list(stages_range):
    stages_dict = {1: 'prepare_data', 2: 'compute & aggregate', 3: 'haplotype analysis & graphs'}
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


def runner(input_dir, reference_file, output_dir, stages_range, max_basecall_iterations, part_size):
    if not max_basecall_iterations:
        max_basecall_iterations = 2   #TODO: what should this be?!
    if not part_size:
        part_size = 10000
    stages = get_stages_list(stages_range)
    log = pipeline_logger(logger_name='AccuNGS-Runner', log_folder=output_dir)
    log.info(f"Running stages: {stages}")
    if 'prepare_data' in stages:
        data_dir = os.path.join(output_dir, "data")
        os.makedirs(data_dir, exist_ok=True)
        log.info("Preparing data")
        prepare_data(input_dir=input_dir, output_dir=data_dir, part_size=part_size, opposing_strings=None)
    else:
        data_dir = input_dir
    if 'compute & aggregate' in stages:
        processing_dir = os.path.join(output_dir, "processing")
        os.makedirs(processing_dir, exist_ok=True)
        data_files = get_files_in_dir(data_dir)
        fastq_files = [file_path for file_path in data_files if "fastq.part_" in os.path.basename(file_path)]
        log.info(f"Running compute & aggregate on {len(fastq_files)} files with params: ")  # TODO: add blast and basecall params
        consensus_file = os.path.join(output_dir, 'consensus_with_indels.fasta')  # TODO: add option to use without indels
        for basecall_iteration_counter in range(1, max_basecall_iterations + 1):
            #TODO: test that it sometimes takes more than just 2 iterations.
            log.info(f"compute & aggregate {basecall_iteration_counter}/{max_basecall_iterations}")
            parallel_compute(processing_dir, fastq_files, reference_file)
            aggregate_computation_output(input_dir=processing_dir, output_dir=output_dir,
                                         reference=reference_file)
            consensus = get_sequence_from_fasta(consensus_file)
            alignment_score = get_alignment_score(consensus=consensus, reference_file=reference_file)
            log.info(f"Iteration done with alignment score: {round(alignment_score,4)}")
            if alignment_score == 1:
                break
            reference_file = consensus_file
        log.info(f"Computation Done. Output files are in directory: {output_dir}")
    if 'haplotype analysis & graphs' in stages:
        freqs_file_path = os.path.join(output_dir, 'freqs.tsv')
        linked_mutations_path = os.path.join(output_dir, 'linked_mutations.tsv')
        mutation_read_list_path = os.path.join(output_dir, 'mutation_read_list.tsv')
        log.info(f"Calculating linked mutations with params: ")  # TODO: add params.
        calculate_linked_mutations(freqs_file_path=freqs_file_path, mutation_read_list_path=mutation_read_list_path,
                                   output=linked_mutations_path)
        stretches = os.path.join(output_dir, 'streches.tsv')
        log.info(f"Aggregating linked mutations to stretches with params: ")
        calculate_stretches(linked_mutations_path, max_pval=10 ** -9, distance=10, output=stretches)  #TODO: refactor that function
        graph_haplotypes(input_file=stretches, number_of_stretches=3, output_dir=output_dir)
        blast_file = os.path.join(output_dir, 'blast.tsv')
        read_counter_file = os.path.join(output_dir, 'read_counter.tsv')
        summary_graphs = os.path.join(output_dir, 'summary.png')
        log.info("Drawing summary graphs!")
        graph_summary(freqs_file=freqs_file_path, blast_file=blast_file, read_counter_file=read_counter_file,
                      stretches_file=stretches, output_file=summary_graphs, min_coverage=3)  # TODO: default min coverage?
    #TODO: test everything, finish up, pbs_runner,

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
    args = parser.parse_args()
    runner(input_dir=args.input_dir, output_dir=args.output_dir, reference_file=args.reference_file,
           stages_range=args.stages_range, max_basecall_iterations=args.max_basecall_iterations,
           part_size=args.part_size)
