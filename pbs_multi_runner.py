import argparse
import json
from pbs_runner import pbs_runner
from utils import get_config


def multi_runner(params_list):
    """
    This runs pbs_runner multiple times and gets its parameters from a list containing the relevant parameters.
    The parameters need to be organised in a list of dictionaries where each dictionary contains the parameters of a
    different run. For example:
    [ {param1: run1_value1, param2: run1_value2 ... },
      {param1: run2_value2, param2: run2_value2, ...},
      ... ]
    If a parameter is not specified it will be retrieved from the defaults described in the config.ini file located in the
    installation directory.

    It is also possible to run this from the cli and supplying the path to a json file containing the said list of
    running parameters.
    """

    for params_dict in params_list:
        args = dict(get_config()['runner_defaults'])                                    # get runner defaults
        args.update(
            {key: value for key, value in dict(get_config()['pbs_defaults']).items()})  # overide with pbs defaults
        args.update(
            {key: value for key, value in params_dict.items() if value is not None})    # overide with params dict
        pbs_runner(input_dir=args['input_dir'], output_dir=args['output_dir'], reference_file=args['reference_file'],
                   max_basecall_iterations=args['max_basecall_iterations'],
                   quality_threshold=args['quality_threshold'], task=args['blast_task'], cleanup=args['cleanup'],
                   evalue=args['blast_evalue'], dust=args['blast_dust'], num_alignments=args['blast_num_alignments'],
                   mode=args['blast_mode'], perc_identity=args['blast_perc_identity'], db_comment=args['db_comment'],
                   min_coverage=args['min_coverage'], soft_masking=args['blast_soft_masking'],
                   with_indels=args['with_indels'], gmem=args['gmem'], stretches_pvalue=args['stretches_pvalue'],
                   stretches_distance=args['stretches_distance'], calculate_haplotypes=args['calculate_haplotypes'],
                   stretches_to_plot=args['stretches_to_plot'], max_read_size=args['stretches_max_read_size'],
                   cpu_count=args['cpu_count'], overlapping_reads=args['overlapping_reads'], db_path=args['db_path'],
                   after_jobid=args['after_jobid'], job_suffix=args['job_suffix'], alias=args['alias'],
                   custom_command=args['custom_command'], queue=args['queue'], default_command=args['default_command'],
                   pbs_cmd_path=args['pbs_cmd_path'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--params_file", required=True, help="Path to config file containing run parameters")
    params_file = parser.parse_args().params_file
    with open(params_file) as jsonfile:
        json_string = jsonfile.read()
        params_list = json.loads(json_string)
    multi_runner(params_list)
