Overview
========

AccuNGS is a sequencing approach aimed at accurate deep sequencing of DNA or RNA,
especially with long genomic regions and low input that often requires amplification. 
It is composed of several techniques that eliminate most errors that are 
created during library prep and sequencing. 

When to use AccuNGS?
^^^^^^^^^^^^^^^^^^^^
There are many sequencing techniques out there. What scenarios are suitable
for AccuNGS sequencing? If your experiment meets the following criteria,
AccuNGS sequencing may be good for you!

* If you want to identify alleles in a population of genomes as rare as 1:10,000.

* If you intend to apply deep sequencing (reach coverage x1,000 and deeper).

* If you don't expect too admixed population of many diverged genomes. 
  AccuNGS identifies variants that are relative to a central consensus sequence, 
  which may be absent when many diverged genomes are present.

.. warning:: 
    Since deep sequencing is usually obtained on short genomic regions, the software
    associated with AccuNGS was optimized to such scenario. For example, it works 
    smoothly on full genomes of viruses. The software is **not designed** to run on 
    long genomes, i.e. the entire human genome. 

AccuNGS's sequencing principles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
The main idea behind this approach is to reduce the error sources to a minimum. 
This is achieved by adhering to the following principles:

* High-yield reverse transcription reactions (RT) [if required]
* High-fidelity PCR reactions [if required]
* High-fidelity tagmentation 
* Size selection for insert size that equals read size
* Illumina paired-end sequencing
* Sequencing of a homogeneous control sample for fitting substitution-specific distributions of errors

.. figure:: _static/Fig1.png
    :scale: 40%
    :align: center
    :alt: AccuNGS principles
    :figclass: align-center
	
    AccuNGS principles.

Exact equipment and materials used are described in :ref:`cite`. 

AccuNGS - code
^^^^^^^^^^^^^^
The code accompanying AccuNGS is divided into several parts, each carrying a different stage
in the computational analysis of sequencing output. The different flows are:

* :doc:`basecall`
* :doc:`variants`
* :doc:`haplotypes`
* :doc:`barcodes`


