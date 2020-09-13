"""
This script takes an input folder containing fastq files, merges them using the new pipeline and than runs both
the perl (old) and python (new) pipeline.
In the output folder there will be the outputs of both pipelines as well as a log file under .log.
Under the 'analysis' folder a joined dataset of both output data will be created as well as several visualisations
showing some main differences between the outputs.
"""

import argparse
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pbs_runner import create_pbs_cmd_file, submit_cmdfile_to_pbs, pbs_runner
from plotting import set_plots_size_params
from utils import get_files_by_extension
STERNLAB_PATH = "/sternadi/home/volume2/ita/AccuNGS-private"


def wrangle_freqs_df(data):
    data = pd.DataFrame.groupby(data, level=[0, 1, 2]).sum()
    data["frequency"] = round(data["base_counter"] / data["coverage"], 6)
    data["probability"] = round(1 - np.power(10, np.log10(1 - data["frequency"] + 1e-07) * (data["coverage"] + 1)), 2)
    # Get ranks
    pd.DataFrame.reset_index(data, level=[0, 1, 2], inplace=True)
    pd.DataFrame.sort_values(data, by=['ref_position', 'frequency', 'base'], ascending=[True, False, False],
                             inplace=True)
    data["coverage_to_set_rank"] = np.where(data["coverage"] > 0, 1, 0)
    data["rank"] = (pd.Series.cumsum(pd.Series(data["coverage_to_set_rank"])) - 1) % 5
    del data["coverage_to_set_rank"]
    # Create freqs file
    pd.DataFrame.set_index(data, keys=["ref_position", "base"], inplace=True)
    return data


def _create_python_output_folder(output_folder):
    python_output_path = os.path.join(output_folder, 'python_output')
    if not os.path.exists(python_output_path):
        os.mkdir(python_output_path)
    return python_output_path


def _create_perl_output_folder(output_folder):
    perl_output_path = os.path.join(output_folder, 'perl_output')
    if not os.path.exists(perl_output_path):
        os.mkdir(perl_output_path)
        perl_tmp_folder = os.path.join(perl_output_path, 'tmp')
        if not os.path.exists(perl_tmp_folder):
            os.mkdir(perl_tmp_folder)  # seems like this is required by the pipeline..
    return perl_output_path


def _get_python_runner_flags(output_folder):
    ret = {}
    ret['i'] = os.path.join(output_folder, 'data')
    ret['o'] = _create_python_output_folder(output_folder)
    return ret


def create_perl_runner_cmdfile(input_data_folder, output_folder, reference_file, alias, pipeline_arguments,
                               python_runner_flags):
    perl_output_path = _create_perl_output_folder(output_folder)
    perl_runner_path = os.path.join(STERNLAB_PATH, 'pipeline_runner.py')
    perl_runner_cmd = f"python {perl_runner_path} -i {python_runner_flags['o']} -o {perl_output_path} " \
                      f"-r {reference_file} -NGS_or_Cirseq 1 -rep {pipeline_arguments['repeats']} " \
                      f"-ev {pipeline_arguments['evalue']} -b {pipeline_arguments['blast']} " \
                      f"-q {pipeline_arguments['q_score']}"
    if len(get_files_by_extension(input_data_folder, 'gz')) > 0:
        perl_runner_cmd = perl_runner_cmd + " -t z"
    """
    the input for perl_runner_cmd is the output of python_runner_cmd because the python pipeline first created
    the fastq files which both pipelines use.
    """
    cmd_file_path = os.path.join(output_folder, 'compare_pipelines.cmd')
    create_pbs_cmd_file(path=cmd_file_path, alias=alias, cmd=perl_runner_cmd, output_logs_dir=output_folder)
    return cmd_file_path


def get_single_freq_file_path(path, freq_file_suffix):
    freq_files = [f for f in os.listdir(path) if f.find(freq_file_suffix) != -1]
    if len(freq_files) == 0:
        raise Exception(f"Could not find file containing {freq_file_suffix} file in {path} !")
    elif len(freq_files) > 1:
        raise Exception(f"Found more than one file containing '{freq_file_suffix}' in {path}..!")
    return os.path.join(path, freq_files[0])


def get_python_freqs(python_output_path):
    freq_file_path = get_single_freq_file_path(python_output_path, "merge.freqs.csv")
    return pd.read_csv(freq_file_path).set_index(['ref_position', 'base'], drop=True)


def get_perl_freqs(perl_output_path):
    """
    Wrangles perl freqs file to fit with Python freqs file.
    """
    perl_freqs_path = get_single_freq_file_path(perl_output_path, '.freqs')
    pe_df = pd.read_csv(perl_freqs_path, sep='\t', index_col=[0, 1, 2], usecols=[0, 1, 2, 3, 4])
    pe_df.reset_index(inplace=True)
    pe_df.columns = ['ref_position', 'base', 'freq', 'ref_base', 'coverage']
    pe_df.set_index(['ref_position', 'ref_base', 'base'], drop=True, inplace=True)
    pe_df['base_counter'] = pe_df['coverage'] * pe_df['freq']
    pe_df.drop(['freq'], axis=1, inplace=True)
    pe_df = wrangle_freqs_df(pe_df)
    return pe_df


def get_freqs_data(output_folder):
    python_output_path = _create_python_output_folder(output_folder)
    py_df = get_python_freqs(python_output_path)
    perl_output_path = _create_perl_output_folder(output_folder)
    pe_df = get_perl_freqs(perl_output_path)
    data = {'py': py_df, 'pe': pe_df}
    return data


def create_analyze_data_cmdfile(output_folder, alias, previous_jobid):
    cmd_file_path = os.path.join(output_folder, 'analyze_data.cmd')
    this_module = os.path.basename(os.path.normpath(os.path.abspath(__file__)))[:-3]
    output_folder_string = '"' + output_folder + '"'
    cmd = f"cd {os.path.join(STERNLAB_PATH, 'Python_pipeline')}; python -c " \
          f"'from {this_module} import analyze_data; analyze_data({output_folder_string})'"  # <- dirty hack for PBS
    create_pbs_cmd_file(cmd=cmd, alias=alias, path=cmd_file_path, output_logs_dir=output_folder,
                        run_after_job_id=previous_jobid)
    return cmd_file_path


def _apply_invert_deletions(row, col):
    if row.base == '-':
        return row[col] * -1
    else:
        return row[col]


def plot_indels(data, output_folder):
    df = data.copy()
    plt.figure(figsize=(20, 10))
    df['base_counter_pe'] = df.apply(lambda row: _apply_invert_deletions(row, 'base_counter_pe'), axis=1)
    df['base_counter_py'] = df.apply(lambda row: _apply_invert_deletions(row, 'base_counter_py'), axis=1)
    indels_pe = df[(df.ref_base_pe == '-') | (df.base == '-')]
    indels_py = df[(df.ref_base_py == '-') | (df.base == '-')]
    plt.figure(figsize=(20, 10))
    plt.plot(indels_pe.index, indels_pe.base_counter_pe, label=f'perl indels', alpha=0.5)
    plt.plot(indels_py.index, indels_py.base_counter_py, label=f'python indels', alpha=0.5)
    plt.xlabel('ref_position')
    plt.ylabel(' deletions <-- base counter --> insertions')
    plt.title('Indels Coverage: perl vs python (positive is insertion, negative is deletion)')
    plt.legend()
    plt.savefig(os.path.join(output_folder, 'indels.png'))


def drop_indels(df):
    return df[(df.base != '-') & (df.ref_base_pe != '-') & (df.ref_base_py != '-')].copy()


def plot_coverage_diff(data, output_folder):
    noindels = drop_indels(data)
    plt.figure(figsize=(20, 10))
    noindels.fillna(0, inplace=True)
    noindels['cov_diff'] = noindels.coverage_pe - noindels.coverage_py
    plt.figure(figsize=(20, 10))
    plt.plot(noindels.index, noindels['cov_diff'])
    plt.xlabel('ref_position')
    plt.ylabel('coverage difference')
    plt.title('Coverage Difference - perl minus python (excluding indels)')
    plt.savefig(os.path.join(output_folder, 'coverage_diff.png'))


def plot_mutations(data, output_folder):
    noindels = drop_indels(data)
    mutations = noindels[(noindels['rank_pe'] > 0) | (noindels['rank_py'] > 0)]
    for base in ['A', 'C', 'G', 'T']:
        mutated_bases = mutations[
            (mutations.base == base) & ((mutations.frequency_pe > 0) | (mutations.frequency_py > 0))]
        plt.figure(figsize=(20, 10))
        plt.scatter(mutated_bases.index, mutated_bases.frequency_pe, alpha=0.5, label='perl pipeline')
        plt.scatter(mutated_bases.index, mutated_bases.frequency_py, alpha=0.5, label='python pipeline')
        plt.xlabel('ref_position')
        plt.ylabel('frequency')
        plt.title(f'X > {base} Mutation Frequency')
        plt.legend()
        plt.savefig(os.path.join(output_folder, f'mutations_{base}.png'))


def plot_frequency_comparison(data, output_folder):
    noindels = drop_indels(data)
    noindels = noindels[(noindels.coverage_py > 5) | (noindels.coverage_pe > 5)]
    noindels.fillna(0, inplace=True)
    plt.figure(figsize=(20, 10))
    plt.scatter(noindels.frequency_pe, noindels.frequency_py)
    plt.xlabel('frequency in perl pipeline')
    plt.ylabel('frequency in python pipeline')
    plt.title("Frequency Python vs Perl where coverage > 5")
    plt.savefig(os.path.join(output_folder, 'frequency_comparison.png'))


def analyze_data(output_folder):
    """
    This function is called by PBS via the cmdfile created by create_analyze_data_cmdfile
    """
    analysis_folder = os.path.join(output_folder, 'analysis')
    if not os.path.isdir(analysis_folder):
        os.mkdir(analysis_folder)
    data = get_freqs_data(output_folder)
    df = data['pe'].join(data['py'], rsuffix='_py', lsuffix='_pe', how='outer').reset_index().set_index('ref_position')
    set_plots_size_params(20)
    plot_indels(data=df, output_folder=analysis_folder)
    plot_coverage_diff(data=df, output_folder=analysis_folder)
    plot_mutations(data=df, output_folder=analysis_folder)
    plot_frequency_comparison(data=df, output_folder=analysis_folder)
    df.to_csv(os.path.join(analysis_folder, 'data.csv'))


def merge_fastq_files(input_data_folder, output_folder, reference_file, pipeline_arguments):
    """ Merge fastq files using the python_runner """
    pbs_runner(input_dir=input_data_folder, output_dir=output_folder,
               reference_file=reference_file, mode="RefToSeq", evalue=pipeline_arguments['evalue'],
               quality_threshold=pipeline_arguments['q_score'], stages_range=1,
               perc_identity=pipeline_arguments['blast'], max_basecall_iterations=1, dust="no",
               num_alignments=1000000, task="blastn")


def main(args):
    input_data_folder = args.input_data_folder
    output_folder = args.output_folder
    reference_file = args.reference_file
    stages = args.stages
    pipeline_arguments = {'blast': args.blast,
                          'evalue': args.evalue,
                          'repeats': args.repeats,
                          'q_score': args.q_score}
    if 'perl' in stages or 'python' in stages:
        python_runner_flags = _get_python_runner_flags(output_folder=output_folder)
        merge_fastq_files(input_data_folder=input_data_folder, output_folder=output_folder,
                          reference_file=reference_file, pipeline_arguments=pipeline_arguments)
    if 'perl' in stages:
        perl_runner_cmd = create_perl_runner_cmdfile(input_data_folder=python_runner_flags['i'],
                                                     output_folder=output_folder,
                                                     reference_file=reference_file, alias='CmpPLPRL',
                                                     pipeline_arguments=pipeline_arguments,
                                                     python_runner_flags=python_runner_flags)
        perl_job_id = submit_cmdfile_to_pbs(perl_runner_cmd)
    if 'python' in stages:
        python_job_id = pbs_runner(input_dir=python_runner_flags['i'], output_dir=python_runner_flags['o'],
                                   reference_file=reference_file, mode="RefToSeq", evalue=pipeline_arguments['evalue'],
                                   quality_threshold=pipeline_arguments['q_score'], stages_range=[1, 3],
                                   perc_identity=pipeline_arguments['blast'], max_basecall_iterations=1, dust="no",
                                   num_alignments=1000000, task="blastn", alias="CmpPLPY")
    if 'analysis' in stages:
        # note that this will only work if both perl and python output already exist.
        # There is no guarantee that the perl command is done except for the fact that the python is way slower...!
        analyze_cmd_path = create_analyze_data_cmdfile(output_folder, alias='CmpPL-Analyze',
                                                       previous_jobid=f"{python_job_id}:{perl_job_id}")
        submit_cmdfile_to_pbs(analyze_cmd_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_data_folder",
                        help="A folder containing the fastq files to run both pipelines on.",
                        required=True)
    parser.add_argument("-o", "--output_folder",
                        help="Where you want the output files to go",
                        required=True)
    parser.add_argument("-r", "--reference_file",
                        required=True)
    parser.add_argument("-s", "--stages", default=['perl', 'python', 'analysis'],
                        help="A list containing any of ['perl', 'python', 'analysis'] default is all of them.")
    parser.add_argument("-b", "--blast", type=int, help=" percent blast id, default=85", default=85)
    parser.add_argument("-ev", "--evalue", type=float, help="E value for blast, default=1e-7", required=False,
                        default=1e-7)
    parser.add_argument("-x", "--repeats", type=int, help="number of repeats, default=1", required=False, default=1)
    parser.add_argument("-q", "--q_score", type=int, help="Q-score cutoff, default=30", required=False, default=30)
    parser.add_argument("-pr", "--please_remove_double_mapping", default='Y', help='Experimental feature to stop double'
                                                                                   ' mapping. Input can be Y or N, '
                                                                                   'default is Y which means it is off')
    args = parser.parse_args()
    main(args)
