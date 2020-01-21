Overview
========

AccuNGS is a sequencing approach aimed at accurate sequencing of DNA or RNA,
especially with long genomic regions and low input. It is composed of several 
different techniques that eliminate most errors that are created during 
library prep and sequencing process. 

AccuNGS's sequencing principles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

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


