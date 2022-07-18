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
## Running on a PBS cluster
todo: runner, pbs_runner, pbs_multi_runner, project_runner & config.ini

# Overview
The pipeline is divided into 4 parts each having it's own .py file: <br>
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
This is where the heavy duty computation takes place and it is a stepwise 
1. run BLAST and basecall on each of the output files of stage I in parallel.<br>
2. compare the resulting consensus alignment with the reference.
Runs BLAST and then basecall on each of the files created by the previous step in parallel. Loops over every nucleotide and decides whether to filter it out based on the given [parameters](#parameters). After running on all files, creates a freqs file and a consensus alignment. Compares the reference and the consensus scores files, if they are identical continue to next stage. else: rerun this step until they are identical or a maximuam threshold set by --max_basecall_iterations / -m is reached.


## Aggregate
Filename: agrregation.py <br>

## Summarize
Filename: summarize.py <br>

This is a work in progress of a python implementation of the entire pipeline.
In order to install the environment: conda env create --file requirements.txt name_of_new_environment
