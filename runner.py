import argparse
import os
import multiprocessing as mp


from data_preparation import prepare_data
from computation import compute
from aggregation import aggregate_computation_output
from utils import get_files_in_dir, get_sequence_from_fasta

# TODO: relative paths
#       progress bar


def parallel_compute(computation_dir, fastq_files, reference, max_cpus=None):
    if not max_cpus:
        max_cpus = mp.cpu_count()
    pool = mp.Pool(processes=max_cpus)
    results = [pool.apply_async(compute, args=(fastq_file, reference, computation_dir, 30, False))
               for fastq_file in fastq_files]
    output = [res.get() for res in results]
    return output


def get_stages_list(stages_range):
    stages_dict = {1: 'prepare_data', 2: 'compute & aggregate'}
    stages = range(stages_range[0], stages_range[1] + 1)
    return [stages_dict[stage] for stage in stages]


def runner(input_dir, reference_file, output_dir, stages_range, max_basecall_iterations):
    if not max_basecall_iterations:
        max_basecall_iterations = 2   #TODO: what should this be?!
    data_dir = os.path.join(output_dir, "data")
    stages = get_stages_list(stages_range)
    if 'prepare_data' in stages:
        os.makedirs(data_dir, exist_ok=True)
        prepare_data(input_dir=input_dir, output_dir=data_dir)
    computation_dir = os.path.join(output_dir, "computation")
    os.makedirs(computation_dir, exist_ok=True)
    if 'compute & aggregate' in stages:
        data_files = get_files_in_dir(data_dir)
        fastq_files = [file_path for file_path in data_files if "fastq.part_" in os.path.basename(file_path)]
        for basecall_iteration_counter in range(1, max_basecall_iterations + 1):
            results = parallel_compute(computation_dir, fastq_files, reference_file)
            aggregate_computation_output(input_dir=os.path.join(computation_dir, "basecall"), output_dir=output_dir,
                                         reference=reference_file)
            consensus_file = os.path.join(output_dir, 'consensus_with_indels.fasta')
            consensus = get_sequence_from_fasta(consensus_file)
            reference = get_sequence_from_fasta(reference_file)
            if consensus == reference:
                break
            reference_file = consensus_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing basecall files")
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-r", "--reference_file", required=True)
    parser.add_argument("-s", "--stages_range", nargs="+", type=int, help="start and end stages separated by spaces")
    parser.add_argument("-m", "--max_basecall_iterations", type=int,
                        help="number of times to run basecall before giving up equalizing reference with consensus")
    args = parser.parse_args()
    runner(input_dir=args.input_dir, output_dir=args.output_dir, reference_file=args.reference_file,
           stages_range=args.stages_range, max_basecall_iterations=args.max_basecall_iterations)
