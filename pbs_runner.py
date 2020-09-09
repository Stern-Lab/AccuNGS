import datetime
import os

from runner import create_runner_parser


def create_pbs_cmd_file(path, alias, output_logs_dir, cmd, queue="adistzachi@power9", gmem=2, ncpus=1, nodes=1,
                        load_python=True, jnums=None, run_after_job_id=None):
    with open(path, 'w') as o:
        o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
        o.write(f"#PBS -q {queue}\n")
        o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
        o.write("#PBS -N " + alias + "\n")
        o.write(f"#PBS -o {output_logs_dir} \n")
        o.write(f"#PBS -e {output_logs_dir} \n")
        o.write(f"#PBS -l select={nodes}:ncpus={ncpus}:mem={gmem}gb\n")
        if jnums:
            if isinstance(jnums, int):
                o.write(f"#PBS -J 1-{str(jnums)} \n\n")
            else:
                o.write(f"#PBS -J {str(jnums[0])}-{str(jnums[1])} \n\n")
        if run_after_job_id is not None and 'adi' in queue:
            o.write("#PBS -W depend=afterok:" + str(run_after_job_id) + ".power9.tau.ac.il\n\n")
        if run_after_job_id is not None and 'dudu' in queue:
            o.write("#PBS -W depend=afterok:" + str(run_after_job_id) + ".power8.tau.ac.il\n\n")
        if load_python:
            o.write("module load python/python-anaconda3.2019.10\n")
        o.write(cmd)
    o.close()


def submit_cmdfile_to_pbs(cmdfile):
    cmd = "/opt/pbs/bin/qsub " + cmdfile
    result = os.popen(cmd).read()
    return result.split(".")[0]


def runner_cmd(input_dir, output_dir, reference_file, stages_range, max_basecall_iterations, part_size,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, soft_masking, min_coverage,
               consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot,
               max_read_size, base_path):
    runner_path = os.path.join(base_path, 'runner.py')
    if not isinstance(stages_range, int):
        stages_range = f"{stages_range[0]} {stages_range[1]}"
    cmd = f"python {runner_path} -i {input_dir} -o {output_dir} -r {reference_file} -s {stages_range}"
    if max_basecall_iterations:
        cmd += f" -m {max_basecall_iterations}"
    if part_size:
        cmd += f" -p {part_size}"
    if quality_threshold:
        cmd += f" -qt {quality_threshold}"
    if task:
        cmd += f" -bt task"
    if evalue:
        cmd += f" -be {evalue}"
    if dust:
        cmd += f" -bd {dust}"
    if num_alignments:
        cmd += f" -bn {num_alignments}"
    if mode:
        cmd += f" -bm {mode}"
    if perc_identity:
        cmd += f" -bp {perc_identity}"
    if soft_masking:
        cmd += f" -bs {soft_masking}"
    if min_coverage:
        cmd += f" -mc {min_coverage}"
    if consolidate_consensus_with_indels:
        cmd += f" -ccwi {consolidate_consensus_with_indels}"
    if stretches_pvalue:
        cmd += f" -sp {stretches_pvalue}"
    if stretches_distance:
        cmd += f" -sd {stretches_distance}"
    if stretches_to_plot:
        cmd += f" -stp {stretches_to_plot}"
    if max_read_size:
        cmd += f" -smrs {max_read_size}"
    return cmd


def pbs_runner_experimental(input_dir, output_dir, reference_file, stages_range, max_basecall_iterations, part_size,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, soft_masking, min_coverage,
               consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot,
               max_read_size, alias, queue, num_of_nucs):
    # TODO: move defaults to a config file
    serial_job_id = None
    haplo_job_id = None
    graph_job_id = None
    jobs_submitted = []
    base_path = os.path.dirname(os.path.abspath(__file__))
    pbs_logs_dir = os.path.join(output_dir, "pbs_logs")
    os.makedirs(pbs_logs_dir, exist_ok=True)
    if isinstance(stages_range, int):
        stages_range = [stages_range, stages_range]
    if min(stages_range) <= 3:
        serial_cmdfile = os.path.join(pbs_logs_dir, f"AccuNGS_123.cmd")
        stages_123 = (stages_range[0], min(3, stages_range[1]))
        cmd = runner_cmd(input_dir=input_dir, output_dir=output_dir, reference_file=reference_file,
                         stages_range=stages_123, max_basecall_iterations=max_basecall_iterations,
                         part_size=part_size, quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust,
                         num_alignments=num_alignments, mode=mode, perc_identity=perc_identity,
                         soft_masking=soft_masking, min_coverage=min_coverage,
                         consolidate_consensus_with_indels=consolidate_consensus_with_indels,
                         stretches_pvalue=stretches_pvalue, stretches_distance=stretches_distance,
                         stretches_to_plot=stretches_to_plot, max_read_size=max_read_size, base_path=base_path)
        alias = "AccuNGS_123"
        create_pbs_cmd_file(serial_cmdfile, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=100,
                            ncpus=40)
        serial_job_id = submit_cmdfile_to_pbs(serial_cmdfile)
        jobs_submitted.append(serial_job_id)
    if (min(stages_range) <= 4) and (max(stages_range) >= 4):
        freqs_file_path = os.path.join(output_dir, 'freqs.tsv')
        alias = "AccuNGS_haplo"
        compute_haplo_path = os.path.join(pbs_logs_dir, f"AccuNGS_4.cmd")
        linked_mutations_dir = os.path.join(output_dir, 'processing', 'linked_mutations')
        mutations_linking_path = os.path.join(base_path, 'mutations_linking.py')
        cmd = f"python {mutations_linking_path} -x $PBS_ARRAY_INDEX -f {freqs_file_path} " \
              f"-r {output_dir}/mutation_read_list.tsv " \
              f"-o {linked_mutations_dir}/$PBS_ARRAY_INDEX_linked_mutations.tsv"
        if max_read_size:
            cmd += f" -m {max_read_size}"
        haplo_logs = os.path.join(pbs_logs_dir, 'linked_mutations_logs')
        os.makedirs(haplo_logs, exist_ok=True)
        create_pbs_cmd_file(compute_haplo_path, alias, output_logs_dir=haplo_logs, cmd=cmd, queue=queue, gmem=20,
                            ncpus=1, jnums=num_of_nucs, run_after_job_id=serial_job_id)
        haplo_job_id = submit_cmdfile_to_pbs(compute_haplo_path)
        jobs_submitted.append(haplo_job_id)
    if 5 in stages_range:
        alias = "AccuNGS_graph"
        graph_haplo_path = os.path.join(pbs_logs_dir, "AccuNGS_5.cmd")
        cmd = runner_cmd(input_dir=input_dir, output_dir=output_dir, reference_file=reference_file,
                         stages_range=5, max_basecall_iterations=max_basecall_iterations,
                         part_size=part_size, quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust,
                         num_alignments=num_alignments, mode=mode, perc_identity=perc_identity,
                         soft_masking=soft_masking, min_coverage=min_coverage,
                         consolidate_consensus_with_indels=consolidate_consensus_with_indels,
                         stretches_pvalue=stretches_pvalue, stretches_distance=stretches_distance,
                         stretches_to_plot=stretches_to_plot, max_read_size=max_read_size, base_path=base_path)
        create_pbs_cmd_file(graph_haplo_path, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=2,
                            ncpus=1, run_after_job_id=haplo_job_id)
        graph_job_id = submit_cmdfile_to_pbs(graph_haplo_path)
        jobs_submitted.append(graph_job_id)
    print(f"Jobs submitted: {jobs_submitted}")
    print(f"Output files will be in {output_dir} ")
    print(f"cmd files and pbs logs under will be in {pbs_logs_dir} ")
    print(f"runner log file will be in {output_dir}.log ")


def pbs_runner(input_dir, output_dir, reference_file, stages_range, max_basecall_iterations, part_size,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, soft_masking, min_coverage,
               consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot,
               max_read_size, alias, queue):
    # TODO: move defaults to a config file
    # TODO: optimize part size depending on input size and number of CPUs
    base_path = os.path.dirname(os.path.abspath(__file__))
    pbs_logs_dir = os.path.join(output_dir, "pbs_logs")
    os.makedirs(pbs_logs_dir, exist_ok=True)
    if isinstance(stages_range, int):
        stages_range = [stages_range, stages_range]
    current_time = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    cmd_path = os.path.join(pbs_logs_dir, f'AccuNGS_{current_time}.cmd')
    cmd = runner_cmd(input_dir=input_dir, output_dir=output_dir, reference_file=reference_file,
                     stages_range=stages_range, max_basecall_iterations=max_basecall_iterations,
                     part_size=part_size, quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust,
                     num_alignments=num_alignments, mode=mode, perc_identity=perc_identity,
                     soft_masking=soft_masking, min_coverage=min_coverage,
                     consolidate_consensus_with_indels=consolidate_consensus_with_indels,
                     stretches_pvalue=stretches_pvalue, stretches_distance=stretches_distance,
                     stretches_to_plot=stretches_to_plot, max_read_size=max_read_size, base_path=base_path)
    create_pbs_cmd_file(cmd_path, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=100,
                        ncpus=50)  #todo: optimize ncpus and put it as a param
    job_id = submit_cmdfile_to_pbs(cmd_path)
    print(f"Job {job_id} submitted.")
    print(f"Output files will be in {output_dir}")
    print(f"cmd file pbs log will be in {pbs_logs_dir}")
    print(f"runner log file will be in {os.path.join(output_dir, '.log')}")


if __name__ == "__main__":
    parser = create_runner_parser()
    parser.add_argument("-a", "--alias", default="AccuNGS", help="job alias visible in qstat (default: AccuNGS)")
    parser.add_argument("-q", "--queue", default="adistzachi@power9",
                        help="PBS queue to run on (default: adistzachi@power9)")
    args = parser.parse_args()
    pbs_runner(input_dir=args.input_dir, output_dir=args.output_dir, reference_file=args.reference_file,
               stages_range=args.stages_range, max_basecall_iterations=args.max_basecall_iterations,
               part_size=args.part_size, quality_threshold=args.quality_threshold, task=args.blast_task,
               evalue=args.blast_evalue, dust=args.blast_dust, num_alignments=args.blast_num_alignments,
               mode=args.blast_mode, perc_identity=args.blast_perc_identity, soft_masking=args.blast_soft_masking,
               min_coverage=args.min_coverage, consolidate_consensus_with_indels=args.consolidate_consensus_with_indels,
               stretches_pvalue=args.stretches_pvalue, stretches_distance=args.stretches_distance, alias=args.alias,
               stretches_to_plot=args.stretches_to_plot, max_read_size=args.stretches_max_read_size, queue=args.queue)
