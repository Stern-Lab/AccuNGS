.. _Sanger: https://en.wikipedia.org/wiki/Phred_quality_score

Overview
========

AccuNGS is a sequencing approach aimed at accurate deep sequencing of DNA or RNA,
especially for clinical samples with low input that are difficult to sequence. 
It is composed of several stages that eliminate most errors that are created 
during library preparataion and sequencing. 

When to use AccuNGS?
^^^^^^^^^^^^^^^^^^^^
There are many sequencing techniques out there. What scenarios are suitable
for AccuNGS sequencing? If your experiment meets the following criteria,
AccuNGS sequencing may be good for you!

* If you want to identify alleles in a population that are as rare as 1:10,000.

* If you plan to sequence deeply (reach coverage of x1,000 and higher).

* If your population is not very heterogeneous of many diverged genomes. 
  AccuNGS identifies variants that are relative to a central consensus sequence.

.. warning:: 
    Since deep sequencing is usually obtained on short genomic regions, the software
    associated with AccuNGS was optimized to such scenario. For example, it works 
    smoothly on full genomes of viruses. The software is **not designed** to run on 
    long genomes, i.e. the entire human genome. 

AccuNGS's sequencing principles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
The main idea behind this approach is to reduce the errors to a minimum. 
This is achieved by adhering to the following principles:

* High-yield reverse transcription reactions (RT) [if required]
* High-fidelity PCR reactions [if required]
* High-fidelity tagmentation 
* Size selection for insert size that equals read size
* Illumina paired-end sequencing
* Sequencing of a homogeneous control sample for fitting specific distributions of errors

.. figure:: _static/Fig1.png
    :scale: 40%
    :align: center
    :alt: AccuNGS principles
    :figclass: align-center
	
    AccuNGS principles.

AccuNGS was calibrated as described in :ref:`cite`. 

AccuNGS - code
^^^^^^^^^^^^^^
The code accompanying AccuNGS is divided into several stages, each carrying out 
a different stage in the computational analysis of sequencing output. 
The different stages are:

* :doc:`basecall`
* :doc:`variants`
* :doc:`haplotypes`
* :doc:`barcodes`

.. note::
    AccuNGS code assumes its sequencing input quality scores to be compatible 
    with the `Sanger`_ format encoded as ASCII +33 (standard Illumina FASTQ 
    output since 2012).
