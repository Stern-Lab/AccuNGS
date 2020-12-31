"""
This is where the magic happens!
This tool ... TODO: what is this and what does it do?!

It is meant to be able to run locally and in pbs_runner.py there is also specific support for pbs cluster systems.

___Overview___
The pipeline is divided into 4 parts each having it's own .py file.
I   -  Preperation (data_preperation.py)
II  -  Processing (processing.py)
III -  Aggregation (aggregation.py)
IV  -  Haplotype Inference (mutations_linking.py)

I - Data Prepeation
Input: directory containing fastq/gz files or a directory containing such directories.
Output: fastq files in sizes ready for efficient processing.

Depending on currently availble RAM and CPUs or given values in -mm / --max_memory and -cc / --cpu_count the script
will divide the files into an efficient number of files to later run the processing script on.
If given an merge_opposing it will also merge the forward and backward reads in corresponding fastq files.

II - Processing
Input: fastq files from part I and a reference fasta file.
Output: called_bases - every base called and its attributes.
        ignored_bases - every base ignored, its attributes and why it was ignored.
        suspicious_reads - every suspicious read and why it is suspicious.
        ignored_reads - every read ignored and why.
        read_counter - every read and a number describing how many times it was aligned by blast.
        blast files - alignment output files.

This is where the main logic of the pipeline happens and the runner runs this on each of the fastq files in parallel.
Basically it first aligns the reads to the reference with blast and then goes over every read and nucleotide and decides
whether to filter them out or leave them in.

III - Aggregation
Input: directory containing output of basecalling of stage II
Output: concatenations - concatenated outputs of stage II
        read_id_prefix_file - a json file containing a dictionary of read prefixes to save memory in the other files.
        mutation_read_list - a file describing in which reads each mutation appeared.
        freqs - a frequencies file describing the different alleles and their frequencies for each position.
        consensus_with_indels - a fasta file of the majority frequency derived from the freqs file including indels
        consensus_without_indels - a fasta file of the majority frequency derived from the freqs file exclusing indels
        read_counter - a file counting how many alignments were called for each read.

IV - Haplotype Inference
Input: directory containing output of stage III
Output: linked_mutations - file containing pairs of mutatations, their pvalues and their frequencies.
        stretches - file containing the aggregations of the stretches by their frequencies.
"""
import argparse
import concurrent.futures
import getpass
import json
import os
import multiprocessing as mp
import shutil

from datetime import datetime
from functools import partial
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
    get_files_by_extension, concatenate_files_by_extension, get_config, md5_dir, md5_file
from haplotypes.co_occurs_to_stretches import calculate_stretches


def parallel_process(processing_dir, fastq_files, reference_file, quality_threshold, task, evalue, dust, num_alignments,
                     soft_masking, perc_identity, mode, reads_overlap):
    process_fastq_partial = partial(process_fastq, reference=reference_file, output_dir=processing_dir,
                                    quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust,
                                    num_alignments=num_alignments, soft_masking=soft_masking,
                                    perc_identity=perc_identity, mode=mode, reads_overlap=reads_overlap)
    try:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executor.map(process_fastq_partial, fastq_files)
    except:
        raise Exception("Processing failed. This may be due to RAM overload."
                        "To avoid this try running with --max_memory max_available_RAM_in_MB")

    path = os.path.join(processing_dir, 'exceptions_during_process_fastq.txt')
    if os.path.exists(path):
        with  open(path, 'r') as file_exception:
            content = file_exception.read()
            file_exception.close()
            exception_printed = "There were exception in the parallel processing. see " \
                            "processing\Exceptions_during_run_blast.txt. "
            if content.find("blastn -out") > 0:
                exception_printed += " There might be a problem involving blast installation"
            if len(content) > 0:
                raise Exception(exception_printed)


def set_filenames(output_dir):
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
    filenames['basecall_dir'] = os.path.join(filenames['processing_dir'], 'basecall')
    os.makedirs(filenames['data_dir'], exist_ok=True)
    os.makedirs(filenames["basecall_dir"], exist_ok=True)
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
    # TODO: optimize memory usage and better exceptions.
    mutation_read_list = pd.read_csv(mutation_read_list_path, sep="\t")
    mutation_read_list['ref_pos'] = mutation_read_list['ref_pos'].round(3)
    mutation_read_list = mutation_read_list.set_index(['ref_pos', 'read_base'], drop=True)
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
        end_position = positions[end_index - 1] + 0.999  # include insertions
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


def get_consensus_path(basecall_iteration_counter, with_indels, iteration_data_dir):
    consensus = os.path.join(iteration_data_dir, "consensus_with")
    if with_indels != "Y":
        consensus += "out"
    consensus += f"_indels_{basecall_iteration_counter}.fasta"
    return consensus


def update_meta_data(output_dir, status, db_path, params=None):
    # TODO: write running time in meta_data
    json_file = os.path.join(output_dir, 'meta_data.json')
    if params is not None:
        meta_data = params
        del meta_data['filenames']
        del meta_data['log']
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
    # TODO: get user feedback on what db functionality should be
    os.makedirs(db_path, exist_ok=True)
    outputs = [f.path for f in os.scandir(db_path) if f.is_dir()]
    db_rows = []
    for directory in outputs:
        json_file = os.path.join(directory, "meta_data.json")
        if os.path.isfile(json_file):
            with open(json_file) as read_handle:
                meta_data = json.load(read_handle)
            db_rows.append(meta_data)
    db = pd.DataFrame.from_dict(db_rows)
    db.to_csv(os.path.join(db_path, 'db.tsv'), sep='\t')


def create_dir_also_if_exists(output_dir):
    try:
        os.makedirs(output_dir, exist_ok=False)
    except FileExistsError:
        i = 2
        while i < 1000:
            try:
                output_dir = f"{output_dir}_{i}"
                os.makedirs(output_dir, exist_ok=False)
            except FileExistsError:
                i += 1
    return output_dir


def assign_output_dir(db_path, alias=None):
    now = datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
    output_dir_name = now
    if alias:
        output_dir_name = alias + "_" + output_dir_name
    output_dir = os.path.join(db_path, output_dir_name)
    output_dir = create_dir_also_if_exists(output_dir)
    return output_dir


def infer_haplotypes(cpu_count, filenames, linked_mutations_dir, log, max_read_size, output_dir, basecall_dir,
                     stretches_distance, stretches_pvalue):
    os.makedirs(linked_mutations_dir, exist_ok=True)
    # TODO: optimize part size
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


def process_data(with_indels, dust, evalue, fastq_files, log, max_basecall_iterations,
                 min_coverage, mode, num_alignments, opposing_strings, output_dir, perc_identity, processing_dir,
                 quality_threshold, reference_file, soft_masking, task, basecall_dir):
    reads_overlap = True if opposing_strings == 'Y' or opposing_strings == 'y' else False
    for basecall_iteration_counter in range(1, max_basecall_iterations + 1):
        log.info(f"Processing fastq files iteration {basecall_iteration_counter}/{max_basecall_iterations}")
        parallel_process(processing_dir=processing_dir, fastq_files=fastq_files, reference_file=reference_file,
                         quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust, mode=mode,
                         num_alignments=num_alignments, soft_masking=soft_masking, perc_identity=perc_identity,
                         reads_overlap=reads_overlap)
        iteration_data_dir = os.path.join(output_dir, 'iteration_data')
        os.makedirs(iteration_data_dir, exist_ok=True)
        alignment_score = check_consensus_alignment_with_ref(reference_file=reference_file,
                                                             iteration_counter=basecall_iteration_counter,
                                                             basecall_dir=basecall_dir,
                                                             with_indels=with_indels,
                                                             iteration_data_dir=iteration_data_dir,
                                                             min_coverage=min_coverage)
        log.info(f'Iteration alignment score: {round(alignment_score, 4)}')
        if alignment_score == 1:
            break
        consensus_path = get_consensus_path(basecall_iteration_counter=basecall_iteration_counter,
                                            with_indels=with_indels, iteration_data_dir=iteration_data_dir)
        reference_file = consensus_path
    return reference_file


def runner(input_dir, reference_file, output_dir, max_basecall_iterations, min_coverage, db_comment,
           quality_threshold, task, evalue, dust, num_alignments, soft_masking, perc_identity, mode, max_read_size,
           with_indels, stretches_pvalue, stretches_distance, stretches_to_plot, cleanup,
           cpu_count, overlapping_reads, db_path, max_memory, calculate_haplotypes="Y"):
    if not db_path:
        db_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'db')
    if not output_dir:
        output_dir = assign_output_dir(db_path)
    log = pipeline_logger(logger_name='AccuNGS-Runner', log_folder=output_dir)
    try:
        filenames = set_filenames(output_dir=output_dir)
        if not cpu_count:
            cpu_count = mp.cpu_count()
        input_dir_hash = md5_dir(input_dir)
        reference_file_hash = md5_file(reference_file)
        params = locals().copy()
        update_meta_data(params=params, output_dir=output_dir, status='Setting up...', db_path=db_path)
        log.debug(f"runner params: {params}")  # TODO: why does this contain status..?
        log.info("Preparing data")
        update_meta_data(output_dir=output_dir, status='Preparing data...', db_path=db_path)
        prepare_data(input_dir=input_dir, output_dir=filenames['data_dir'], overlapping_reads=overlapping_reads,
                     cpu_count=cpu_count, max_memory=max_memory)
        data_files = get_files_in_dir(filenames['data_dir'])
        fastq_files = [file_path for file_path in data_files if "fastq.part_" in os.path.basename(file_path)]
        log.info(f"Processing {len(fastq_files)} fastq files.")
        update_meta_data(output_dir=output_dir, status='Processing data...', db_path=db_path)
        reference_file = process_data(with_indels=with_indels, dust=dust,
                                      evalue=evalue, fastq_files=fastq_files, log=log, soft_masking=soft_masking,
                                      max_basecall_iterations=max_basecall_iterations, min_coverage=min_coverage,
                                      mode=mode, num_alignments=num_alignments, opposing_strings=overlapping_reads,
                                      output_dir=output_dir, perc_identity=perc_identity, reference_file=reference_file,
                                      processing_dir=filenames['processing_dir'], quality_threshold=quality_threshold,
                                      task=task, basecall_dir=filenames['basecall_dir'])
        log.info("Aggregating processed fastq files outputs...")
        aggregate_processed_output(input_dir=filenames['processing_dir'], output_dir=output_dir,
                                   reference=reference_file, min_coverage=min_coverage)
        log.info("Generating graphs...")
        graph_summary(freqs_file=filenames['freqs_file_path'], blast_file=filenames['blast_file'],
                      read_counter_file=filenames['read_counter_file'], stretches_file=filenames['stretches'],
                      output_file=filenames['summary_graphs'], min_coverage=min_coverage,
                      stretches_to_plot=stretches_to_plot)  # TODO: drop low quality mutations?
        log.info(f"Most outputs are ready in {output_dir} !")
        if calculate_haplotypes == "Y" or calculate_haplotypes == "y":
            log.info(f"Calculating linked mutations...")
            update_meta_data(output_dir=output_dir, status='Inferring haplotypes...', db_path=db_path)
            infer_haplotypes(cpu_count=cpu_count, filenames=filenames,
                             linked_mutations_dir=filenames['linked_mutations_dir'],
                             log=log, max_read_size=max_read_size, output_dir=output_dir,
                             stretches_pvalue=stretches_pvalue,
                             basecall_dir=filenames['basecall_dir'], stretches_distance=stretches_distance)
            graph_summary(freqs_file=filenames['freqs_file_path'], blast_file=filenames['blast_file'],
                          read_counter_file=filenames['read_counter_file'], stretches_file=filenames['stretches'],
                          output_file=filenames['summary_graphs'], min_coverage=min_coverage,
                          stretches_to_plot=stretches_to_plot)  # TODO: drop low quality mutations?
            graph_haplotypes(input_file=filenames['stretches'], number_of_stretches=stretches_to_plot,
                             output_dir=output_dir)
        if cleanup == "Y":
            log.info(f"Removing directory {filenames['basecall_dir']} to save space.")
            shutil.rmtree(filenames['basecall_dir'])
            with open(filenames['basecall_dir'], 'w') as f:
                f.write("Directory 'basecall' which contains the original parts of the output files has been removed. "
                        "To prevent it from being remove run with option --cleanup N")
        update_meta_data(output_dir=output_dir, status='Done', db_path=db_path)
        log.info(f"Done!")
    except Exception as e:
        log.exception(e)
        update_meta_data(output_dir=output_dir, status="Failed! see logs for details.", db_path=db_path)


def create_runner_parser():
    # TODO: dynamic defaults?
    parser = argparse.ArgumentParser(description="Note: Default parameter values are retrieved from config.ini in your "
                                                 "installation directory.")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing fastq/gz files or sub directories containg fastq/gz files.")
    parser.add_argument("-o", "--output_dir", help="A directory for output files. "
                                                   "If none is given will put it in the db")
    parser.add_argument("-r", "--reference_file", required=True, help="Reference file to align against.")
    parser.add_argument("-m", "--max_basecall_iterations", type=int,
                        help="Number of times to rerun with previous consensus as the new reference before giving up.")
    parser.add_argument("-or", "--overlapping_reads",
                        help="Y/N, merge opposing reads in the same directory. This assumes 2 fastq/gz files in each "
                             "sub directory of the input_dir and would drop all non overlapping areas of the reads.")
    parser.add_argument("-bt", "--blast_task", help="blast's task parameter")
    parser.add_argument("-be", "--blast_evalue", help="blast's evalue parameter", type=float)
    parser.add_argument("-bd", "--blast_dust", help="blast's dust parameter")
    parser.add_argument("-bn", "--blast_num_alignments", type=int, help="blast's num_alignments parameter")
    parser.add_argument("-bp", "--blast_perc_identity", type=int, help="blast's perc_identity parameter")
    parser.add_argument("-bs", "--blast_soft_masking", help="blast's soft_masking parameter")
    parser.add_argument("-bm", "--blast_mode", help="RefToSeq or SeqToRef")  # TODO: docs
    parser.add_argument("-qt", "--quality_threshold", type=int,
                        help="phred score must be higher than this to be included")
    parser.add_argument("-mc", "--min_coverage", type=int,
                        help="Minimal coverage required for a position to be considered as consensus")
    parser.add_argument("-wi", "--with_indels", help="Y/N, create consensus with or without indels")
    parser.add_argument("-sp", "--stretches_pvalue", type=float,
                        help="only consider joint mutations with pvalue below this value")
    parser.add_argument("-sd", "--stretches_distance", type=float,
                        help="mean transitive distance between joint mutations to calculate stretches")
    parser.add_argument("-stp", "--stretches_to_plot", type=int,
                        help="number of stretches to plot in deep dive")
    parser.add_argument("-smrs", "--stretches_max_read_size", type=int,
                        help="look this many positions forward for joint mutations")
    parser.add_argument("-c", "--cleanup", help="Remove processing/basecall directory when done")
    parser.add_argument("-cc", "--cpu_count", help="max number of cpus to use (None means all)", type=int)
    parser.add_argument("-db", "--db_path", help='path to db directory')
    parser.add_argument("-dbc", "--db_comment", help='comment to store in db')
    parser.add_argument("-mm", "--max_memory", help='limit memory usage to this many megabytes '
                                                    '(None would use available memory when starting to run)')
    parser.add_argument("-ch", "--calculate_haplotypes", help='Y/N, Run pipeline including calculating haplotypes')
    return parser


if __name__ == "__main__":
    parser = create_runner_parser()
    parser_args = vars(parser.parse_args())
    args = dict(get_config()['runner_defaults'])
    args.update({key: value for key, value in parser_args.items() if value is not None})
    runner(input_dir=args['input_dir'], output_dir=args['output_dir'], reference_file=args['reference_file'],
           max_basecall_iterations=int(args['max_basecall_iterations']), overlapping_reads=args['overlapping_reads'],
           quality_threshold=int(args['quality_threshold']), task=args['blast_task'], max_memory=args['max_memory'],
           evalue=float(args['blast_evalue']), dust=args['blast_dust'],
           num_alignments=int(args['blast_num_alignments']),
           mode=args['blast_mode'], perc_identity=float(args['blast_perc_identity']), cpu_count=args['cpu_count'],
           min_coverage=int(args['min_coverage']), db_comment=args['db_comment'],
           soft_masking=args['blast_soft_masking'],
           stretches_pvalue=float(args['stretches_pvalue']), stretches_distance=float(args['stretches_distance']),
           cleanup=args['cleanup'], with_indels=args['with_indels'], calculate_haplotypes=args['calculate_haplotypes'],
           stretches_to_plot=int(args['stretches_to_plot']), max_read_size=int(args['stretches_max_read_size']),
           db_path=args['db_path'])
