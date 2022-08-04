This is a work in progress of a python implementation of the entire pipeline.
In order to install the environment: <br>
<code>conda env create --file requirements.txt name_of_new_environment</code>

# <center>___AccuNGS___
# Introduction
<b><i>AccuNGS</b></i> is the main computational pipeline used by SternLab. <br>
This Python implementation overides and improves over the previous perl implementation which is still available under branch perl.<br>
Within AccuNGS wet protocol, the sequencing library is created to maximize overlap of the Forward and Reverse reads of a paired-end Illumina sequencing run, in a manner that allows for two observations of each base of the original insert for the purpose of increasing the accuracy of basecalling.

# Installation
todo: gotta package and test this thing.. 

# Usage
AccuNGS pipeline was designed to run locally while using available memory and cpus efficiently and it also natively supports running on a [PBS cluster](#running-on-a-pbs-cluster) system. AccuNGS has three [required parameters](#required-parameters) and all other parameters have default values which can be changed by editing the file <i>config.ini</i> in the installation directory. AccuNGS pipeline was designed to run on fastq files with maximal overlap between the forward and reverse reads but can also be run on fastq files without any overlap using the using the parameter <code>--overlapping_reads N</code>

## Output
todo: explain freqs, consensus, called_bases, filttered things and maybe examples of graphs.

## Parameters
#### Required Parameters
  - <code>-i / --input_dir</code> - Path to directory containing fastq/fastq.gz files or sub directories containg fastq/fastq.gz files.
- <code>-o / --output_dir</code> - Path to directory where output files will go.
- <code>-r / --reference_file</code> - Full path to reference fasta file which fastq files will be aligned against using BLAST.

#### Basecall Parameters
- <code>-m / --max_bascall_interations</code> - Maximum number of Process loop iterations (see [Process](#process)).
- <code>-or / --overlapping_reads</code> - Y/N/P, run pipeline with, without or with partial overlapping reads. Y would ignore all bases without an overlap, N assumes reads are independent and ignores overlap completely and P uses all available information.
- <code>-qt / --quality_threshold</code> - Filter out all nucleotides with phred score lower than this.
- <code>-mc / --min_coverage</code> - Positions with less than this coverage will be replaced by N's in the consensus.
- <code>-mf / --min_frequency</code> - Positions with less than this frequency will be replaced by N's in the consensus.
- <code>-ar / --align_to_ref</code> - Y/N, align consensus to original reference.

#### BLAST Parameters
The pipeline allows control of some of BLAST's parameters in order to get a better alignment. For an indepth view of all of BLAST's parameters see [BLAST's documentation](https://www.ncbi.nlm.nih.gov/books/NBK279684/).
- <code>-bt / --blast_task</code> - BLAST's task parameter.
- <code>-be / --blast_evalue</code> - BLAST's evalue parameter.
- <code>-bd / --blast_dust</code> - BLAST's dust parameter.
- <code>-bn / --blast_num_alignments</code> - BLAST's num_alignments parameter.
- <code>-bp / --blast_perc_identity</code> - BLAST's perc_identity parameter.
- <code>-bs / --blast_soft_masking</code> - BLAST's soft_masking parameter.

#### Efficieny Parameters
- <code>-c / --cleanup</code> - Remove redundant basecall directory when done.
- <code>-cc / --cpu_count</code> - Max number of CPUs to use (0 means all).
- <code>-mm / --max_memory</code> - Limit memory usage to this many megabytes (0 would use available memory in the begining of execution).
  

## Running on a PBS cluster
todo: runner, pbs_runner, pbs_multi_runner, project_runner & config.ini

# Algorithm Overview
The pipeline is divided into 4 stages each having it's own .py file: <br>
I   -  [Prepare Data](#prepare-data)  - Prepare files for efficient parallel processing <br>
II  -  [Process](#process)  - Parallel run BLAST and basecall on each of the outputs of the previous step <br>
III -  [Aggregate](#aggregate) - Aggregate outputs of parallel runs (aggregation.py) <br>
IV  -  [Summarize](#summarize) - Draw some graphs and a short text summary <br>

## Prepare Data 
Filename: data_preperation.py <br>
Takes a directory containing fastq/gz files or a directory containing sub directories containing fastq/gz files and outputs fastq files which are ready for processing. <br>
If parameter --overlapping_reads / -or is set to 'Y' or 'P' (partial) the opposing reads will be merged into a single file. <br>
Based on the values in --max_memory / --mm and --cpu_count / -cc the pipeline will divivde the input files into an efficient number of files to run simultanesouly. If the values are not given the pipeline will try to estimate these numbers according to the currently available system resources.

## Process 
Filename: processing.py <br>
Runs [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and then basecall on each of the files created in stage I in parallel. BLAST is a local alignment tool used to align the reads with the reference file given by parameter --reference / -r . All the relevant BLAST variables can be set with their [corresponding parameters](#blast-parameters). Basecall loops over every nucleotide and decides whether to filter it out based on the relevant [parameters](#basecall-parameters). After the parallel run on all stage I output files is complete, a freqs file and a consensus alignment are created. If the concensus and the reference files are identical, the pipeline continues to stage III. If they are not identical, stage II is repeated but with the newly created consensus file given as the new reference file. This loop continues until the consensus converges fully or a maximuam threshold of iterations set by --max_basecall_iterations / -m is reached.

## Aggregate
Filename: agrregation.py <br>
Aggregates and cleans up outputs of stage II and creates the main output files.
  
## Summarize
Filename: summarize.py <br>
Creates a graphical and textual summary of the output.
