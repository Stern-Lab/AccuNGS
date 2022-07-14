"""
Run pbs_multi_runner with identical parameters on multiple datasets arranged in sub-directories of input_dir.
output_dir will contain directories corresponding to a pipeline ran on each sub-directory in input_dir.
"""
import os
import json
from pbs_runner import get_pbs_args
from pbs_multi_runner import multi_runner


def create_params_list(args):
    parent_input = args['input_dir']
    parent_output = args['output_dir']
    params_list = []
    dir_list = [x for x in os.listdir(parent_input) if os.path.isdir(os.path.join(parent_input, x))]
    if len(dir_list) == 0:
        raise Exception(f"input_dir ({parent_input}) does not contain sub-directories."
                        f"Input data should be arranged in sub-directories of input_dir.")
    for dir_name in dir_list:
        dir_path = os.path.join(parent_input, dir_name)
        if os.path.isdir(dir_path):
            params_dict = args.copy()
            params_dict['input_dir'] = dir_path
            params_dict['output_dir'] = os.path.join(parent_output, dir_name)
            params_list.append(params_dict)
    return params_list


def run_project(args):
    output_dir = args['output_dir']
    os.makedirs(output_dir, exist_ok=True)
    params_list_file = os.path.join(output_dir, 'project_params.json')
    params_list = create_params_list(args)
    with open(params_list_file, 'w') as write_handle:
        json.dump(params_list, write_handle, indent=4)
    print(f"Project params file in {params_list_file}")
    job_ids = multi_runner(params_list)
    return job_ids


if __name__ == "__main__":
    args = get_pbs_args()
    run_project(args)
