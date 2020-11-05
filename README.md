This is the repository accompanying the AccuNGS paper (PLoS Pathog 2020; DOI 10.1371/journal.ppat.1009029) and is subject to the AccuNGS paper user license. 

Within AccuNGS protocol, the sequencing library is created to maximize overlap of the Forward and Reverse reads of a paired-end Illumina sequencing run, in a manner that allows for two observations of each base of the original insert. This practically cancels much of the sequencing by synthesis errors. In order to support paired-end-aware analysis, we have created a basecalling pipeline. 

The usage of the this pipeline should usually follow this pattern:
1. Concatenate all reads in the Forward and Reverse files, with multiple N bases as separators (script: merger.py)
2. Map reads to reference genome by BLAST command line. 
3. Basecall bases with sufficient quality score (script: base_call_and_freqs.pl)

Parts 2-3 can be run in parallel on a cluster of computational nodes. We provide a similar pipeline compatible with PBS. 

If a homogeneous control was sequenced in parallel to the sample(s), it is useful to identify variants by using AccuNGS variant caller, which assigns p-values to the observed frequencies based on the probability they represent true variation. The variant caller outputs a file similar to the output provided by V-phaser2 (script: AccuNGS_variant_caller.py). (V-phaser2 website: https://www.broadinstitute.org/viral-genomics/v-phaser-2 )

AccuNGS also supports degenerate barcodes recovery and analysis (also known as primer-ID, if barcode is added during RT of RNA). Code for barcode extraction and analysis is provided (folder: barcode_analysis)

We also provide calculators to accompany the AccuNGS protocol, for sequencing experiment planning.
