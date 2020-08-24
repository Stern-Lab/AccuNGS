import argparse
import os
import multiprocessing as mp
from data_preparation import prepare_data
from computation import compute
from aggregation import aggregate_computation_output
from utils import get_files_by_extension, get_files_in_dir

# TODO: relative paths
#       progress bar

def mp_compute(computation_dir, fastq_files, reference, max_cpus=None):
    if not max_cpus:
        max_cpus = mp.cpu_count()
    pool = mp.Pool(processes=max_cpus)
    results = [pool.apply_async(compute, args=(fastq_file, reference, computation_dir, 30, False))
               for fastq_file in fastq_files]
    output = [res.get() for res in results]
    return output


def runner(input_dir, reference, output_dir, prepare_data_stage=False, computation_stage=False, aggregation_stage=True,
           consolidate_consensus=True, max_basecall_iterations=2):
    data_dir = os.path.join(output_dir, "data")

    if prepare_data_stage:
        os.makedirs(data_dir, exist_ok=True)
        prepare_data(input_dir=input_dir, output_dir=data_dir)
    computation_dir = os.path.join(output_dir, "computation")
    os.makedirs(computation_dir, exist_ok=True)
    if consolidate_consensus:
        run_1 = os.path.join(computation_dir, "run_1")
        os.makedirs(run_1, exist_ok=True)
        computation_dir = run_1
    data_files = get_files_in_dir(data_dir)
    fastq_files = [file_path for file_path in data_files if "fastq.part_" in os.path.basename(file_path)]
    results = mp_compute(computation_dir, fastq_files, reference)
    if aggregation_stage:
        aggregate_computation_output(input_dir=os.path.join(computation_dir, "basecall"), output_dir=output_dir)
    if consolidate_consensus:
        for basecall_iteration_counter in range(2, max_basecall_iterations + 1):
            new_ref = make_reference_from_freqs(freqs, reference)
            if reference == new_ref:
                break
            iteration_output_dir = os.path.join(computation_dir, f"run_{basecall_iteration_counter}")
            results = mp_compute(iteration_output_dir, fastq_files, new_ref)
            reference = new_ref
            new_ref = make_reference_from_freqs(os.path.join(iteration_output_dir, "freqs.tsv"))
        if reference != new_ref:
            log.warning(f"Could not consolidate consensus in {max_basecall_iterations} iterations!")






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Path to directory containing basecall files")
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-r", "--reference", required=True)
    args = parser.parse_args()
    runner(input_dir=args.input_dir, output_dir=args.output_dir, reference=args.reference)
