import argparse
import concurrent.futures
import getpass
import json
import os
import multiprocessing as mp
import shutil

from datetime import datetime
from functools import partial
from coolname import generate_slug
import pandas as pd
from Bio import pairwise2

from data_preparation import prepare_data
from graph_haplotypes import graph_haplotypes
from mutations_linking import get_variants_list, get_mutations_linked_with_position
from plotting import graph_summary
from processing import process_fastq
from aggregation import aggregate_processed_output, create_freqs_file, create_mutation_read_list_file
from logger import pipeline_logger
from utils import get_files_in_dir, get_sequence_from_fasta, get_mp_results_and_report, create_new_ref_with_freqs, \
    get_files_by_extension, concatenate_files_by_extension, get_config
from haplotypes.co_occurs_to_stretches import calculate_stretches


def parallel_process(processing_dir, fastq_files, reference_file, quality_threshold, task, evalue, dust, num_alignments,
                     soft_masking, perc_identity, mode, reads_overlap):
    process_fastq_partial = partial(process_fastq,  reference=reference_file, output_dir=processing_dir,
                                    quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust,
                                    num_alignments=num_alignments, soft_masking=soft_masking,
                                    perc_identity=perc_identity, mode=mode, reads_overlap=reads_overlap)
    try:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executor.map(process_fastq_partial, fastq_files)
    except:
        raise Exception("Processing failed. This may be due to RAM overload."
                        "To avoid this try running with --max_memory max_available_RAM_in_MB")


def set_filenames(output_dir, db_path):
    if not output_dir:
        output_dir = assign_output_dir(db_path)
    filenames = {"freqs_file_path": os.path.join(output_dir, 'freqs.tsv'),
                 "linked_mutations_path": os.path.join(output_dir, 'linked_mutations.tsv'),
                 "mutation_read_list_path": os.path.join(output_dir, 'mutation_read_list.tsv'),
                 "stretches": os.path.join(output_dir, 'stretches.tsv'),
                 "blast_file": os.path.join(output_dir, 'blast.tsv'),
                 "read_counter_file": os.path.join(output_dir, 'read_counter.tsv'),
                 "summary_graphs": os.path.join(output_dir, 'summary.png'),
                 "processing_dir": os.path.join(output_dir, "processing"),
                 'linked_mutations_dir': os.path.join(output_dir, "linked_mutations"),
                 'data_dir': os.path.join(output_dir, "data")}
    os.makedirs(filenames['data_dir'], exist_ok=True)
    os.makedirs(filenames["processing_dir"], exist_ok=True)
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
    #TODO: optimize memory usage and proper exceptions.
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


def update_meta_data(output_dir, status, db_path, params=None):
    json_file = os.path.join(output_dir, 'meta_data.json')
    if params is not None:
        meta_data = params
        del meta_data['filenames']
        meta_data['username'] = getpass.getuser()
        meta_data['start_time'] = datetime.now().strftime('%Y-%m-%d-%H:%M')
    else:
        with open(json_file) as read_handle:
            meta_data = json.load(read_handle)
    meta_data['status'] = status
    with open(json_file, 'w') as write_handle:
        json.dump(meta_data, write_handle)
    build_db(db_path)


def build_db(db_path):
    outputs = [f.path for f in os.scandir(db_path) if f.is_dir()]
    db_rows = []
    for dir in outputs:
        json_file = os.path.join(dir, "meta_data.json")
        if os.path.isfile(json_file):
            with open(json_file) as read_handle:
                meta_data = json.load(read_handle)
            db_rows.append(meta_data)
    db = pd.DataFrame.from_dict(db_rows)
    db.to_csv(os.path.join(db_path, 'db.tsv'), sep='\t')


def assign_output_dir(db_path):
    now = datetime.now().strftime('%Y-%m-%d-%H:%M')
    random_name = generate_slug(2)
    output_dir_name = random_name + "_" + now
    output_dir = os.path.join(db_path, output_dir_name)
    return output_dir


def infer_haplotypes(cpu_count, filenames, linked_mutations_dir, log, max_read_size, output_dir, processing_dir,
                     stretches_distance, stretches_pvalue):
    os.makedirs(linked_mutations_dir, exist_ok=True)
    # TODO: optimize part size
    basecall_dir = os.path.join(processing_dir, 'basecall')
    called_bases_files = get_files_by_extension(basecall_dir, "called_bases")
    mutation_read_list_path = os.path.join(output_dir, "mutation_read_list.tsv")
    create_mutation_read_list_file(called_bases_files=called_bases_files, output_path=mutation_read_list_path)
    parallel_calc_linked_mutations(freqs_file_path=filenames['freqs_file_path'], cpu_count=cpu_count,
                                   mutation_read_list_path=filenames['mutation_read_list_path'],
                                   output_dir=linked_mutations_dir, max_read_length=max_read_size,
                                   part_size=100)  # TODO: drop low quality mutations?, set part_size as param.
    log.info(f"Aggregating linked mutations to stretches...")
    concatenate_files_by_extension(input_dir=linked_mutations_dir, extension='tsv',
                                   output_path=filenames['linked_mutations_path'])
    calculate_stretches(filenames['linked_mutations_path'], max_pval=stretches_pvalue, distance=stretches_distance,
                        output=filenames['stretches'])  # TODO: refactor that function


def process_data(consolidate_consensus_with_indels, dust, evalue, fastq_files, log, max_basecall_iterations,
                 min_coverage, mode, num_alignments, opposing_strings, output_dir, perc_identity, processing_dir,
                 quality_threshold, reference_file, soft_masking, task):
    reads_overlap = bool(opposing_strings)
    for basecall_iteration_counter in range(1, max_basecall_iterations + 1):
        log.info(f"Processing fastq files iteration {basecall_iteration_counter}/{max_basecall_iterations}")
        # TODO: whats up with the different modes?!?
        parallel_process(processing_dir=processing_dir, fastq_files=fastq_files, reference_file=reference_file,
                         quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust, mode=mode,
                         num_alignments=num_alignments, soft_masking=soft_masking, perc_identity=perc_identity,
                         reads_overlap=reads_overlap)
        iteration_data_dir = os.path.join(output_dir, 'iteration_data')
        os.makedirs(iteration_data_dir, exist_ok=True)
        alignment_score = check_consensus_alignment_with_ref(reference_file=reference_file,
                                                             iteration_counter=basecall_iteration_counter,
                                                             basecall_dir=os.path.join(processing_dir, 'basecall'),
                                                             with_indels=consolidate_consensus_with_indels,
                                                             iteration_data_dir=iteration_data_dir,
                                                             min_coverage=min_coverage)
        log.info(f'Iteration alignment score: {round(alignment_score, 4)}')
        if alignment_score == 1:
            break
        consensus_path = get_consensus_path(basecall_iteration_counter, consolidate_consensus_with_indels,
                                            iteration_data_dir)
        reference_file = consensus_path
    return reference_file


def runner(input_dir, reference_file, output_dir, max_basecall_iterations, min_coverage, db_comment,
           quality_threshold, task, evalue, dust, num_alignments, soft_masking, perc_identity, mode, max_read_size,
           consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot, cleanup,
           cpu_count, opposing_strings, db_path, max_memory):
    try:
        filenames = set_filenames(output_dir=output_dir, db_path=db_path)
        if not cpu_count:
            cpu_count = mp.cpu_count()
        params = locals().copy()
        update_meta_data(params=params, output_dir=output_dir, status='Setting up...', db_path=db_path)
        log = pipeline_logger(logger_name='AccuNGS-Runner', log_folder=output_dir)
        log.debug(f"runner params: {params}")
        log.info("Preparing data")
        update_meta_data(output_dir=output_dir, status='Preparing data...', db_path=db_path)
        prepare_data(input_dir=input_dir, output_dir=filenames['data_dir'], overlap_notation=opposing_strings,
                     cpu_count=cpu_count, max_memory=max_memory)
        data_files = get_files_in_dir(filenames['data_dir'])
        fastq_files = [file_path for file_path in data_files if "fastq.part_" in os.path.basename(file_path)]
        log.info(f"Processing {len(fastq_files)} fastq files.")
        update_meta_data(output_dir=output_dir, status='Processing data...', db_path=db_path)
        reference_file = process_data(consolidate_consensus_with_indels=consolidate_consensus_with_indels, dust=dust,
                                      evalue=evalue, fastq_files=fastq_files, log=log, soft_masking=soft_masking,
                                      max_basecall_iterations=max_basecall_iterations, min_coverage=min_coverage,
                                      mode=mode, num_alignments=num_alignments, opposing_strings=opposing_strings,
                                      output_dir=output_dir, perc_identity=perc_identity, reference_file=reference_file,
                                      processing_dir=filenames['processing_dir'], quality_threshold=quality_threshold,
                                      task=task)
        log.info("Aggregating processed fastq files outputs...")
        aggregate_processed_output(input_dir=filenames['processing_dir'], output_dir=output_dir, cleanup=cleanup,
                                   reference=reference_file, min_coverage=min_coverage)
        log.info("Generating graphs...")
        graph_summary(freqs_file=filenames['freqs_file_path'], blast_file=filenames['blast_file'],
                      read_counter_file=filenames['read_counter_file'], stretches_file=filenames['stretches'],
                      output_file=filenames['summary_graphs'], min_coverage=min_coverage,
                      stretches_to_plot=stretches_to_plot)  # TODO: drop low quality mutations?
        log.info(f"Most outputs are ready in {output_dir} !")
        log.info(f"Calculating linked mutations...")
        update_meta_data(output_dir=output_dir, status='Inferring haplotypes...', db_path=db_path)
        infer_haplotypes(cpu_count=cpu_count, filenames=filenames, linked_mutations_dir=filenames['linked_mutations_dir'],
                         log=log, max_read_size=max_read_size, output_dir=output_dir, stretches_pvalue=stretches_pvalue,
                         processing_dir=filenames['processing_dir'], stretches_distance=stretches_distance)
        graph_summary(freqs_file=filenames['freqs_file_path'], blast_file=filenames['blast_file'],
                      read_counter_file=filenames['read_counter_file'], stretches_file=filenames['stretches'],
                      output_file=filenames['summary_graphs'], min_coverage=min_coverage,
                      stretches_to_plot=stretches_to_plot)  # TODO: drop low quality mutations?
        graph_haplotypes(input_file=filenames['stretches'], number_of_stretches=stretches_to_plot,
                         output_dir=output_dir)
        if cleanup == "Y":
            shutil.rmtree(filenames['processing_dir'])
        update_meta_data(output_dir=output_dir, status='Done', db_path=db_path)
        log.info(f"Done!")
    except Exception as e:
        log.exception(e)
        update_meta_data(output_dir=output_dir, status="Failed! see logs for details.", db_path=db_path)

    #TODO: test everything, finish up, drop cpu_count?




def create_runner_parser():
    # TODO: present dynamic defaults from config file
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing basecall files")
    parser.add_argument("-o", "--output_dir")
    parser.add_argument("-r", "--reference_file", required=True)
    parser.add_argument("-m", "--max_basecall_iterations", type=int, default=1,
                        help="number of times to run basecall before giving up equalizing reference with consensus")
    parser.add_argument("-on", "--overlap_notation", nargs="+", type=str,
                        help="Notation of overlapping reads in the same directory to merge (default: '_R1 _R2') and "
                             "passing nothing would run without overlap")
    parser.add_argument("-bt", "--blast_task", help="blast's task parameter (default: blastn")
    parser.add_argument("-be", "--blast_evalue", help="blast's evalue parameter (default: 1e-7)", type=float)
    parser.add_argument("-bd", "--blast_dust", help="blast's dust parameter (default: on)")
    parser.add_argument("-bn", "--blast_num_alignments", type=int,
                        help="blast's num_alignments parameter (default: 1000000)")
    parser.add_argument("-bp", "--blast_perc_identity", type=int, help="blast's perc_identity parameter (default: 85)")
    parser.add_argument("-bs", "--blast_soft_masking", help="blast's soft_masking parameter (default: F)")
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
    parser.add_argument("-db", "--db_path", help='path to db directory')
    parser.add_argument("-dbc", "--db_comment", help='comment to store in db')
    parser.add_argument("-mm", "--max_memory", help='limit memory usage to this many megabytes (default: None)')
    return parser


if __name__ == "__main__":
    parser = create_runner_parser()
    parser_args = vars(parser.parse_args())
    args = dict(get_config()['runner_defaults'])
    args.update({key: value for key, value in parser_args.items() if value is not None})
    runner(input_dir=args['input_dir'], output_dir=args['output_dir'], reference_file=args['reference_file'],
           max_basecall_iterations=int(args['max_basecall_iterations']),
           quality_threshold=int(args['quality_threshold']), task=args['blast_task'], max_memory=args['max_memory'],
           evalue=float(args['blast_evalue']), dust=args['blast_dust'], num_alignments=int(args['blast_num_alignments']),
           mode=args['blast_mode'], perc_identity=float(args['blast_perc_identity']),
           min_coverage=int(args['min_coverage']), db_comment=args['db_comment'], soft_masking=args['blast_soft_masking'],
           stretches_pvalue=float(args['stretches_pvalue']), stretches_distance=float(args['stretches_distance']),
           cleanup=args['cleanup'], consolidate_consensus_with_indels=args['consolidate_consensus_with_indels'],
           stretches_to_plot=int(args['stretches_to_plot']), max_read_size=int(args['stretches_max_read_size']),
           cpu_count=args['cpu_count'], opposing_strings=args['overlap_notation'], db_path=args['db_path'])

