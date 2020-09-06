import argparse
import os

from runner import get_stages_list, create_runner_parser
from utils import create_pbs_cmd_file, submit_cmdfile_to_pbs


def pbs_runner(input_dir, output_dir, reference_file, stages_range, max_basecall_iterations, part_size,
               quality_threshold, task, evalue, dust, num_alignments, mode, perc_identity, soft_masking, min_coverage,
               consolidate_consensus_with_indels, stretches_pvalue, stretches_distance, stretches_to_plot,
               max_read_size, alias, queue):
    # TODO: write cmds
    stages = get_stages_list(stages_range)
    serial_stages = ['prepare data', 'process data', 'aggregate output']
    serial_job_id = None
    haplo_job_id = None
    graph_job_id = None
    pbs_logs_dir = os.path.join(output_dir, "pbs_logs")
    os.makedirs(pbs_logs_dir, exist_ok=True)
    if set(serial_stages) & set(stages):  # TODO: something clearer?
        serial_cmdfile = os.path.join(pbs_logs_dir, f"serial_pipeline_stages.cmd")
        cmd = "runner stages 1-3"
        alias += "_serial"
        create_pbs_cmd_file(serial_cmdfile, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=100,
                            ncpus=30, nodes=1, load_python=True, jnum=False, run_after_job_id=None)
        serial_job_id = submit_cmdfile_to_pbs(serial_cmdfile)
    if "compute haplotypes" in stages:
        # TODO: get jnum and write cmd
        alias += "_haplo"
        compute_haplo_path = os.path.join(pbs_logs_dir, f"compute_haplotypes_array.cmd")
        cmd = "runner stage 4 only as array"
        create_pbs_cmd_file(compute_haplo_path, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=100,
                            ncpus=1, nodes=1, load_python=True, jnum=jnum, run_after_job_id=serial_job_id)
        haplo_job_id = submit_cmdfile_to_pbs(compute_haplo_path)
    if "graph haplotypes" in stages:
        alias += "_graph"
        graph_haplo_path = os.path.join(pbs_logs_dir, "graph_haplotypes.cmd")
        cmd = "runner stage 5 only"
        create_pbs_cmd_file(graph_haplo_path, alias, output_logs_dir=pbs_logs_dir, cmd=cmd, queue=queue, gmem=2, ncpus=1,
                            nodes=1, load_python=True, jnum=False, run_after_job_id=haplo_job_id)
        graph_job_id = submit_cmdfile_to_pbs(graph_haplo_path)
    print(f"Jobs submitted: {serial_job_id, haplo_job_id, graph_job_id}")
    print(f"Output files will be in {output_dir} .")
    print(f"cmd files and pbs logs under will be in {pbs_logs_dir} .")
    print(f"runner log file will be in {output_dir}/.log ")


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
