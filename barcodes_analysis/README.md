The general process for extracting barcodes is as the following:
split (merged) fastq -> run alignment in parallel using a computational cluster-> merge outputs -> analyze barcodes distribution

Linux command for splitting original (merged) fastq file(s) to 50k reads per file:
split -l 200000 --numeric-suffixes=1 --suffix-length=3 ${fastq_filename} ${output_prefix_with_dot_in_the_end}

Linux command for rename split files to remove extra "0" in the filenames (after CD'ing into the directory)
for file in `ls`; do a=`echo ${file} | cut -f2 -d'.'`; b=$((10#$a)); mv $file ${output_prefix_with_dot_in_the_end}${b}; done

run on cluster barcode_aligner.py with the barcode fasta (example code attached)
concatenate all outputs to a single file
run analyze_primer_id.py on the concatenated file
