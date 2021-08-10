import os
import time
from datetime import datetime
import json
from time import sleep

from pbs_multi_runner import multi_runner

# TODO: this should not be part of the final repo


def read_json_file(filepath):
    with open(filepath) as jsonfile:
        json_string = jsonfile.read()
        data = json.loads(json_string)
    return data


def wait_till_timeout(start_time, timeout):
    if time.time() - start_time > timeout:
        raise TimeoutError('took too long for tests to run!! check it yourself, loser11!!!')
    sleep(5)


def tester():
    timeout = 120
    params_file = '/sternadi/home/volume2/ita/sternlab-public/tester/params.json'
    params_list = read_json_file(params_file)
    now = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    output_dir = '/sternadi/home/volume2/ita/sternlab-public/tester/tests/' + now
    os.makedirs(output_dir, exist_ok=True)
    for param_dict in params_list:
        param_dict['output_dir'] = os.path.join(output_dir, param_dict['output_dir'])
    multi_runner(params_list)
    print('Giving the jobs a bit of time to start...')
    meta_data_files = [os.path.join(param_dict['output_dir'], 'meta_data.json') for param_dict in params_list]
    start_time = time.time()
    while not os.path.isfile(meta_data_files[0]):
        wait_till_timeout(start_time, timeout)
    print('Finally wrote some metadata, yay!')
    for meta_data_file in meta_data_files:
        meta_data = read_json_file(meta_data_file)
        start_time = time.time()
        while '...' in meta_data['status']:
            wait_till_timeout(start_time, timeout)
        status = meta_data['status'][:6]
        job_name = os.path.dirname(meta_data['output_dir'])
        if 'Done' in status:
            print(f"job {job_name} Done!")
        else:
            print(f'error in job {job_name} !')


if __name__ == "__main__":
    tester()
