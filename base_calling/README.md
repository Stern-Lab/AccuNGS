The general flow of the base calling process in AccuNGS is as the following (LINUX env required):

1. Forward (R1) and Reverse (R2) of a paired-end (Illumina) sequencing run are concatenated using "merger.py" python script. The script simply concatenates the two reads to a single read, with multiple "N" bases between the two.

2. BLAST (v2.7.1+) the merged fastq against the reference. Parameters:

	ref_genome - reference genome (FASTA format)
	
	in_fastq - merged fastq file
	
	max_num_alignments - maximum number of alignments (typically 10x of the length of the input file)
	
	pcID_blast - percents identity of each alignment to the reference. suggested default: 85
	
	out_blast - blast output file
	
	Typical use case:
	
	makeblastdb -in ${in_fastq} -dbtype nucl
	
	blastn -query ${ref_genome} -task blastn -db ${in_fastq} -outfmt \"6 sseqid qstart qend qstrand sstart send sstrand length btop\" -num_alignments ${max_num_alignments} -dust no -soft_masking F -perc_identity ${pcID_blast} -evalue 1e-10 -out ${out_blast};
	
3. Basecalling with "base_call_and_freqs.pl" perl script. Parameters:
	out_blast - result output of the previous BLAST run
	in_fastq - merged fastq file 
	ref_genome - the reference FASTA used for BLAST.
	output_file - the output file name (better use .FREQS for convenience)
	do_gaps - Y if the base calling should report indels  
	min_qual_score - minimum average quality score (on the two reads) to be reported. Typically 30 or 38   
	 
	Typical use case:
	perl base_call_and_freqs.pl ${out_blast} ${in_fastq} ${ref_genome} ${output_file} ${do_gaps} ${min_qual_score}
	      
For large FASTQ files, it is extremely useful to split the initial input FASTQ files into smaller parts and running all code in parallel on a computational grid, and then performing a smart merge on the result FREQS files. A PBS-compatible pipeline is available under 'PBS' directory ("runner.pl"). 
