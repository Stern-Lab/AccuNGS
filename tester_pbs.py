import os
from datetime import datetime
import json
from time import sleep

from pbs_multi_runner import multi_runner

# TODO: this should not be part of the final repo


def tester():
    params_file = '/sternadi/home/volume2/ita/sternlab-public/tester_params.json'
    with open(params_file) as jsonfile:
        json_string = jsonfile.read()
        params_list = json.loads(json_string)
    now = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    output_dir = '/sternadi/home/volume2/ita/sternlab-public/db/tester/' + now
    os.makedirs(output_dir, exist_ok=True)
    for param_dict in params_list:
        param_dict['output_dir'] = os.path.join(output_dir, param_dict['output_dir'])
    multi_runner(params_list)
    print('Waiting for jobs to finish...')
    sleep(10)
    print(f"found {len(os.listdir())} files!")


if __name__ == "__main__":
    tester()
