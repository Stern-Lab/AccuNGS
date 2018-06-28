#!/bin/bash
#PBS -S /bin/bash
#PBS -j oe
#PBS -r y
#PBS -q my_queue
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH 
#PBS -N find_primer_ids
#PBS -l mem=4000mb
#PBS -J 1-2

python barcode_aligner.py barcode.fasta reads_file.$PBS_ARRAY_INDEX output_file_prefix.$PBS_ARRAY_INDEX.barcodes.txt
