import argparse
import json
from pbs_runner import pbs_runner
from utils import get_config

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--params_file", required=True, help="Path to config file containing run parameters")
    params_file = parser.parse_args().params_file
    with open(params_file) as jsonfile:
        json_string = jsonfile.read()
        params_dict = json.loads(json_string)
    for alias, parameters in params_dict.items():
        args = dict(get_config()['runner_defaults'])  # get runner defaults
        args.update({key: value for key, value in dict(get_config()['pbs_defaults']).items()})  # overide with pbs defaults
        args.update({key: value for key, value in parameters.items() if value is not None})  # overide with params file args
        pbs_runner(input_dir=args['input_dir'], output_dir=args['output_dir'], reference_file=args['reference_file'],
                   max_basecall_iterations=args['max_basecall_iterations'],
                   quality_threshold=args['quality_threshold'], task=args['blast_task'], cleanup=args['cleanup'],
                   evalue=args['blast_evalue'], dust=args['blast_dust'], num_alignments=args['blast_num_alignments'],
                   mode=args['blast_mode'], perc_identity=args['blast_perc_identity'], db_comment=args['db_comment'],
                   min_coverage=args['min_coverage'], soft_masking=args['blast_soft_masking'],
                   consolidate_consensus_with_indels=args['consolidate_consensus_with_indels'],
                   stretches_pvalue=args['stretches_pvalue'], stretches_distance=args['stretches_distance'],
                   stretches_to_plot=args['stretches_to_plot'], max_read_size=args['stretches_max_read_size'],
                   cpu_count=args['cpu_count'], overlap_notation=args['overlap_notation'], db_path=args['db_path'],
                   after_jobid=args['after_jobid'], job_suffix=args['job_suffix'], alias=alias,
                   custom_command=args['custom_command'], queue=args['queue'], default_command=args['default_command'])
