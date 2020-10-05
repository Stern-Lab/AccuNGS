import os
from random import randint

from runner import create_runner_parser
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
            o.write(default_command)
        if custom_command:
            o.write(custom_command)
        o.write(cmd)
    o.close()


def submit_cmdfile_to_pbs(cmdfile):
    cmd = "/opt/pbs/bin/qsub " + cmdfile
    result = os.popen(cmd).read()
    return result.split(".")[0]


def runner_cmd(input_dir, output_dir, reference_file, stages_range, max_basecall_iterations, db_path,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, soft_masking, min_coverage,
               consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot,
               max_read_size, base_path, cleanup, cpu_count, overlap_notation):
    runner_path = os.path.join(base_path, 'runner.py')
    if not isinstance(stages_range, int):
        stages_range = f"{stages_range[0]} {stages_range[1]}"
    cmd = f"python {runner_path} -i {input_dir} -o {output_dir} -r {reference_file} -s {stages_range}"
    if max_basecall_iterations:
        cmd += f" -m {max_basecall_iterations}"
    if quality_threshold:
        cmd += f" -qt {quality_threshold}"
    if task:
        cmd += f" -bt {task}"
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
    if cleanup:
        cmd += f" -c {cleanup}"
    if cpu_count:
        cmd += f" -cc {cpu_count}"
    if overlap_notation:
        cmd += f" -on {overlap_notation}"
    if db_path:
        cmd += f" -db {db_path}"
    return cmd


def pbs_runner(input_dir, output_dir, reference_file, stages_range, max_basecall_iterations,db_path,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, overlap_notation,
               soft_masking, min_coverage, consolidate_consensus_with_indels, stretches_pvalue, stretches_distance,
               stretches_to_plot, max_read_size, alias, queue, cleanup, cpu_count, custom_command=None, after_jobid=None,
               job_suffix=None):
    config = get_config()
    defaults = config['pbs_defaults']
    base_path = os.path.dirname(os.path.abspath(__file__))
    pbs_logs_dir = os.path.join(output_dir, "pbs_logs")
    os.makedirs(pbs_logs_dir, exist_ok=True)
    if isinstance(stages_range, int):
        stages_range = [stages_range, stages_range]
    cmd_identifier = randint(42, 777)  # so that we can easily connect cmdfile and job
    alias += f"_{cmd_identifier}"
    cmd_path = os.path.join(pbs_logs_dir, f'{alias}.cmd')
    cmd = runner_cmd(input_dir=input_dir, output_dir=output_dir, reference_file=reference_file,
                     stages_range=stages_range, max_basecall_iterations=max_basecall_iterations,
                     quality_threshold=quality_threshold, task=task, evalue=evalue, dust=dust,
                     num_alignments=num_alignments, mode=mode, perc_identity=perc_identity,
                     soft_masking=soft_masking, min_coverage=min_coverage, cleanup=cleanup, cpu_count=cpu_count,
                     consolidate_consensus_with_indels=consolidate_consensus_with_indels, db_path=db_path,
                     stretches_pvalue=stretches_pvalue, stretches_distance=stretches_distance,
                     stretches_to_plot=stretches_to_plot, max_read_size=max_read_size, base_path=base_path,
                     overlap_notation=overlap_notation)
    default_command = defaults.get('default_command', fallback=None)
    create_pbs_cmd_file(cmd_path, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=100,
                        ncpus=cpu_count, run_after_job_id=after_jobid, job_suffix=job_suffix, custom_command=custom_command,
                        default_command=default_command)  # todo: optimize cpu_count
    job_id = submit_cmdfile_to_pbs(cmd_path)
    print(f"Submitted jod '{alias}' with id {job_id}")
    print(f"Output files will be in {output_dir}")
    print(f"cmd file pbs log will be in {pbs_logs_dir}")
    print(f"runner log file will be in {os.path.join(output_dir, '.log')}")
    return job_id


if __name__ == "__main__":
    parser = create_runner_parser()
    # TODO: pbs-runner defaults, run with config file
    parser.add_argument("-a", "--alias", default="AccuNGS", help="job alias visible in qstat (default: AccuNGS)")
    parser.add_argument("-q", "--queue", default="adistzachi@power9",
                        help="PBS queue to run on (default: adistzachi@power9)")
    parser.add_argument("-j", "--after_jobid",
                        help="Run after successfully completing this jobid")
    parser_args = vars(parser.parse_args())
    args = dict(get_config()['runner_defaults'])
    args.update({key: value for key, value in dict(get_config()['pbs_defaults']).items() if value is not None})
    args.update({key: value for key, value in parser_args.items() if value is not None})
    pbs_runner(input_dir=args['input_dir'], output_dir=args['output_dir'], reference_file=args['reference_file'],
               stages_range=args['stages_range'], max_basecall_iterations=args['max_basecall_iterations'],
               quality_threshold=args['quality_threshold'], task=args['blast_task'],
               evalue=args['blast_evalue'], dust=args['blast_dust'], num_alignments=args['blast_num_alignments'],
               mode=args['blast_mode'], perc_identity=args['blast_perc_identity'], soft_masking=args['blast_soft_masking'],
               min_coverage=args['min_coverage'], consolidate_consensus_with_indels=args['consolidate_consensus_with_indels'],
               stretches_pvalue=args['stretches_pvalue'], stretches_distance=args['stretches_distance'], cleanup=args['cleanup'],
               stretches_to_plot=args['stretches_to_plot'], max_read_size=args['stretches_max_read_size'],
               cpu_count=args['cpu_count'], overlap_notation=args['overlap_notation'], db_path=args['db_path'],
               after_jobid=args['after_jobid'], job_suffix=args['job_suffix'], alias=args['alias'],
               custom_command=args['custom_command'], queue=args['queue'])
