"""
This script runs the runner on a PBS cluster system by creating a pbs cmd file and running it with qsub.
In the output directory there will be an extra directory called pbs_logs which will contain the cmd file and the
output logs of the cluster.
You can set defaults for the pbs_runner in the config.ini file in the installation directory.
These can also be parameters for the runner itself which will be used only when running the pbs_runner.
"""

import os
from random import randint

from runner import create_runner_parser, assign_output_dir
from utils import get_config


def create_pbs_cmd_file(path, alias, output_logs_dir, cmd, queue, gmem=10, ncpus=50, nodes=1, custom_command=None,
                        jnums=None, run_after_job_id=None, job_suffix=None, default_command=None):
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
        if run_after_job_id:
            if job_suffix:
                run_after_job_id = str(run_after_job_id) + job_suffix
            o.write("#PBS -W depend=afterok:" + str(run_after_job_id) + "\n")
        if default_command:
            o.write(default_command + "\n")
        if custom_command:
            o.write(custom_command + "\n")
        o.write(cmd)
    o.close()


def submit_cmdfile_to_pbs(cmdfile, pbs_cmd_path):
    cmd = f"{pbs_cmd_path} {cmdfile}"
    result = os.popen(cmd).read()
    return result.split(".")[0]


def runner_cmd(input_dir, output_dir, reference_file, max_basecall_iterations, db_path, db_comment,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, soft_masking, min_coverage,
               align_to_ref, stretches_pvalue, stretches_distance, stretches_to_plot, python_path, min_frequency,
               max_read_size, base_path, cleanup, cpu_count, overlapping_reads, calculate_haplotypes):
    runner_path = os.path.join(base_path, 'runner.py')
    cmd = f"{python_path} {runner_path} -i {input_dir} -o {output_dir} -r {reference_file} "
    if max_basecall_iterations is not None:
        cmd += f" -m {max_basecall_iterations}"
    if quality_threshold is not None:
        cmd += f" -qt {quality_threshold}"
    if task is not None:
        cmd += f" -bt {task}"
    if evalue is not None:
        cmd += f" -be {evalue}"
    if dust is not None:
        cmd += f" -bd {dust}"
    if num_alignments is not None:
        cmd += f" -bn {num_alignments}"
    if mode is not None:
        cmd += f" -bm {mode}"
    if perc_identity is not None:
        cmd += f" -bp {perc_identity}"
    if soft_masking is not None:
        cmd += f" -bs {soft_masking}"
    if min_coverage is not None:
        cmd += f" -mc {min_coverage}"
    if min_frequency is not None:
        cmd += f" -mf {min_frequency}"
    if align_to_ref is not None:
        cmd += f" -ar {align_to_ref}"
    if stretches_pvalue is not None:
        cmd += f" -sp {stretches_pvalue}"
    if stretches_distance is not None:
        cmd += f" -sd {stretches_distance}"
    if stretches_to_plot is not None:
        cmd += f" -stp {stretches_to_plot}"
    if max_read_size is not None:
        cmd += f" -smrs {max_read_size}"
    if cleanup is not None:
        cmd += f" -c {cleanup}"
    if cpu_count is not None:
        cmd += f" -cc {cpu_count}"
    if overlapping_reads is not None:
        cmd += f" -or {overlapping_reads}"
    if db_path is not None:
        cmd += f" -db {db_path}"
    if db_comment is not None:
        cmd += f" -dbc '{db_comment}'"
    if calculate_haplotypes is not None:
        cmd += f" -ch '{calculate_haplotypes}'"
    return cmd


def pbs_runner(input_dir, output_dir, reference_file, max_basecall_iterations, db_path, db_comment, pbs_cmd_path,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, overlapping_reads, gmem,
               soft_masking, min_coverage, align_to_ref, stretches_pvalue, stretches_distance, calculate_haplotypes,
               stretches_to_plot, max_read_size, alias, queue, cleanup, cpu_count, python_path, min_frequency,
               custom_command=None, after_jobid=None, job_suffix=None, default_command=None):
    # TODO: reintroduce stages into the runner,
    #       do relevant aggregation (when needed) in aggregation instead of linked mutations,
    #       create a dir with files representing mutation_read_list_parts,
    #       use it in a pbs array calling linked mutations on each part
    #       where each job uses 24 cpus or less and limit to 500 jobs per array.

    if not output_dir:
        output_dir = assign_output_dir(db_path, alias)
    base_path = os.path.dirname(os.path.abspath(__file__))
    pbs_logs_dir = os.path.join(output_dir, "pbs_logs")
    os.makedirs(pbs_logs_dir, exist_ok=True)
    if not alias:
        cmd_identifier = randint(42, 777)  # so that we can easily connect cmdfile and job
        alias = f"AccuNGS_{cmd_identifier}"
    cmd_path = os.path.join(pbs_logs_dir, f'{alias}.cmd')
    cmd = runner_cmd(input_dir=input_dir, output_dir=output_dir, reference_file=reference_file,
                     max_basecall_iterations=max_basecall_iterations, max_read_size=max_read_size,
                     quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust, python_path=python_path,
                     num_alignments=num_alignments, mode=mode, perc_identity=perc_identity, db_comment=db_comment,
                     soft_masking=soft_masking, min_coverage=min_coverage, cleanup=cleanup, cpu_count=cpu_count,
                     align_to_ref=align_to_ref, db_path=db_path, calculate_haplotypes=calculate_haplotypes,
                     stretches_pvalue=stretches_pvalue, stretches_distance=stretches_distance, base_path=base_path,
                     stretches_to_plot=stretches_to_plot, overlapping_reads=overlapping_reads,
                     min_frequency=min_frequency)
    create_pbs_cmd_file(cmd_path, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=gmem,
                        ncpus=cpu_count, run_after_job_id=after_jobid, job_suffix=job_suffix,
                        custom_command=custom_command, default_command=default_command)
    job_id = submit_cmdfile_to_pbs(cmd_path, pbs_cmd_path)
    if job_id:
        print(f"Submitted jod '{alias}' with id {job_id}")
        print(f"Output files will be in {output_dir}")
        print(f"runner log file will be in {os.path.join(output_dir, '.log')}")
    else:
        print(f"Could not submit job {alias} to queue!")
    print(f"cmd file and pbs logs are in {pbs_logs_dir}")
    return job_id


def get_pbs_args():
    parser = create_runner_parser()
    parser.add_argument("-a", "--alias", help="job alias visible in qstat")
    parser.add_argument("-q", "--queue", help="PBS queue to run on")
    parser.add_argument("-j", "--after_jobid", help="Run after successfully completing this jobid")
    parser.add_argument("-gm", "--gmem", help="Memory in GB to ask for in cmd file")
    parser_args = vars(parser.parse_args())
    args = dict(get_config()['runner_defaults'])  # get runner defaults
    args.update({key: value for key, value in dict(get_config()['pbs_defaults']).items()})  # overide with pbs defaults
    args.update({key: value for key, value in parser_args.items() if value is not None})  # overide with cli args
    return args


if __name__ == "__main__":
    args = get_pbs_args()
    pbs_runner(input_dir=args['input_dir'], output_dir=args['output_dir'], reference_file=args['reference_file'],
               max_basecall_iterations=args['max_basecall_iterations'], custom_command=args['custom_command'],
               quality_threshold=args['quality_threshold'], task=args['blast_task'], db_comment=args['db_comment'],
               evalue=args['blast_evalue'], dust=args['blast_dust'], num_alignments=args['blast_num_alignments'],
               mode=args['blast_mode'], perc_identity=args['blast_perc_identity'],
               min_coverage=args['min_coverage'], cleanup=args['cleanup'], default_command=args['default_command'],
               align_to_ref=args['align_to_ref'], queue=args['queue'], python_path=args['python_path'],
               stretches_pvalue=args['stretches_pvalue'], stretches_distance=args['stretches_distance'],
               stretches_to_plot=args['stretches_to_plot'], max_read_size=args['stretches_max_read_size'],
               cpu_count=args['cpu_count'], overlapping_reads=args['overlapping_reads'], db_path=args['db_path'],
               after_jobid=args['after_jobid'], job_suffix=args['job_suffix'], alias=args['alias'],
               calculate_haplotypes=args['calculate_haplotypes'], pbs_cmd_path=args['pbs_cmd_path'], gmem=args['gmem'],
               soft_masking=args['blast_soft_masking'], min_frequency=args['min_frequency'])
