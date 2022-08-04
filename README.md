# <center>___AccuNGS___
# Introduction
<b><i>AccuNGS</b></i> is the main computational pipeline used by SternLab. <br>
This Python implementation overides and improves over the previous perl implementation which is still available under branch perl.<br>
Within AccuNGS wet protocol, the sequencing library is created to maximize overlap of the Forward and Reverse reads of a paired-end Illumina sequencing run, in a manner that allows for two observations of each base of the original insert for the purpose of increasing the accuracy of basecalling.

# Installation
todo: gotta package and test this thing.. 

# Usage
## Input
todo: supported input files, what is overlap and what is the reference
## Output
todo: explain freqs, consensus, called_bases, filttered things and maybe examples of graphs.
## Parameters
todo: parse --help in here
#### BLAST Parameters
#### Basecall Parameters
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

## Summarize
Filename: summarize.py <br>

This is a work in progress of a python implementation of the entire pipeline.
In order to install the environment: conda env create --file requirements.txt name_of_new_environment
