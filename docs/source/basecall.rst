.. _Python (3.5+): https://www.python.org/downloads/
.. _Perl (5.26+): https://www.perl.org/get.html
.. _Blast (2.2+): https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _BioPython: https://biopython.org/
.. _GitHub repository: https://github.com/SternLabTAU/AccuNGS/

Base calling 
============
The script files supporting base calling are present in the ``base_calling/`` folder in the AccuNGS `GitHub repository`_. 
The input for this stage is a pair of FASTQ files the correspond to an Illumina paired-end sequencing run, and a reference FASTA file.
The output is a file that contains frequencies of different alleles observed for each loci, including insertions and deletions. 
Unmapped loci (i.e. due to lack of coverage) are omitted from the output file.

.. note:: 
    If you have PBS-compatible computational grid, you may want to use the 
    AccuNGS base calling flow from ``base_calling/PBS/`` folder. It will
    streamline the scripts used for this stage, and will parallel the creation 
    of the output by splitting the input for smaller files. Run (``runner.pl``)
    on the merged fastq folder.


Base calling stage
##################

The base calling process in AccuNGS is composed of the following steps, 
that need to be executed one after another.

Merging paired reads
^^^^^^^^^^^^^^^^^^^^
Match forward (R1) and reverse (R2) reads of an Illumina paired-end sequencing
output to each other using ``base_calling/merger.py`` python script. The script simply
concatenates the two reads to a single read, with multiple "N" bases between
the two. In order to run it requires `Python (3.5+)`_ with `BioPython`_ package installed.

.. code-block:: bash

  python base_calling/merger.py ${forward_fastq} ${reverse_fastq} ${merged_fastq}

Turn the merged FASTQ into FASTA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For the next steps, it is required to generate a FASTA file representing the reads. 
A linux command to create it:

.. code-block:: bash

  cat ${merged_fastq} | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${in_fasta}

Run BLAST against the reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This step performs mapping of the reads against the reference sequence, using `Blast (2.2+)`_. 
This step requires the following parameters:	

===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
ref_genome            Text           reference genome (FASTA format)
--------------------- -------------- --------------------------------
in_fasta              Text           merged reads file
--------------------- -------------- --------------------------------
max_num_alignments    Integer        maximum number of alignments (typically 10x the number of reads in the input file)
--------------------- -------------- --------------------------------
pcID_blast            Float          percents identity of each alignment to the reference. suggested default: 85 for lab-derived sequences of 40 for samples from natural populations
--------------------- -------------- --------------------------------
out_blast             Text           BLAST output file
===================== ============== ================================

A typical use case involves creating a BLAST database and then performing the BLAST search:

.. code-block:: bash
    :dedent: 4

        makeblastdb -in ${in_fasta} -dbtype nucl
    	blastn -query ${ref_genome} -task blastn -db ${in_fasta} -outfmt "6 sseqid qstart qend qstrand sstart send sstrand length btop" -num_alignments ${max_num_alignments} -dust no -soft_masking F -perc_identity ${pcID_blast} -evalue 1e-10 -out ${out_blast}

   
Run base-calling script
^^^^^^^^^^^^^^^^^^^^^^^
Base calling is done using ``base_call_and_freqs.pl`` perl script. It requires `Perl (5.26+)`_ installed. 

The perl script requires the following parameters:

===================== ============== ================================ 
Parameter name        Type           Description
===================== ============== ================================
out_blast             Text           result output of the previous BLAST run
--------------------- -------------- --------------------------------
in_fastq              Text           the merged FASTQ file 
--------------------- -------------- --------------------------------
ref_genome            Text           the reference genome (FASTA format) used for BLAST
--------------------- -------------- --------------------------------
output_file           Text           the output file name (better use .FREQS for convenience)
--------------------- -------------- --------------------------------
do_gaps               Text           Y if the base calling should report indels; N otherwise
--------------------- -------------- --------------------------------
min_qual_score        Float          minimum average quality score (on the two reads) to be reported. Typically 30 or 38
===================== ============== ================================

A typical use case:

.. code-block:: bash

  perl base_calling/base_call_and_freqs.pl ${out_blast} ${in_fastq} ${ref_genome} ${output_file} ${do_gaps} ${min_qual_score}

.. _basecall_output:

Output format
#############
The output of this chain of scripts is similar to this 
:download:`example output file <examples/example_output.freqs>`, 
and contains frequencies of different alleles observed at each locus, 
including insertions and deletions (if ``base_call_and_freqs.pl`` was 
run with ``do_gaps=Y``). 

===================== ================ ================================ 
Column name           Type             Description
===================== ================ ================================
Pos                   Float            The locus which this data refers to, when using the reference for numbering. Insertions are marked by consecutive significant figures for this position, for instance insertion after locus ``15`` will be marked as ``15.1``.
--------------------- ---------------- --------------------------------
Base                  String           The identified base. One of {A,C,G,T,-}
--------------------- ---------------- --------------------------------
Freq                  Float            The frequency that this base was observed at this position
--------------------- ---------------- --------------------------------
Ref                   String           The referece base 
--------------------- ---------------- --------------------------------
Coverage              Integer          The number of reads supporting this locus. The number of reads supporting this base are ``Coverage x Freq``
--------------------- ---------------- --------------------------------
Rank                  Integer          The rank of ``Base`` starting at 0 for most prevalent base and descending
===================== ================ ================================

Example: sample data
####################
Here are two example files for 
:download:`Forward (R1) <examples/example_input_S1_L001_R1_001.fastq.gz>` and 
:download:`Reverse (R2) <examples/example_input_S1_L001_R2_001.fastq.gz>` FASTQ files.
Together with a :download:`Reference FASTA file <examples/example_reference.fasta>` 
the base calling flow can be executed, to output the following 
:download:`output file <examples/example_output.freqs>`.

.. code-block:: bash

  gunzip example_input_S1_L001_R1_001.fastq.gz 
  gunzip example_input_S1_L001_R2_001.fastq.gz
  
.. code-block:: bash

  python base_calling/merger.py example_input_S1_L001_R1_001.fastq example_input_S1_L001_R2_001.fastq example_input_S1.merged.fastq

.. code-block:: bash

  cat example_input_S1.merged.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > merged.fasta

.. code-block:: bash
    :dedent: 4

        makeblastdb -in merged.fasta -dbtype nucl
    	blastn -query example_reference.fasta -task blastn -db merged.fasta -outfmt "6 sseqid qstart qend qstrand sstart send sstrand length btop" -num_alignments 100000 -dust no -soft_masking F -perc_identity 85 -evalue 1e-10 -out merged.fasta.out.blast

.. code-block:: bash

  perl base_calling/base_call_and_freqs.pl merged.fasta.out.blast merged.fastq example_reference.fasta merged.fasta.freqs Y 30

The output of this chain of scripts is this :download:`example output file <examples/example_output.freqs>`. 

