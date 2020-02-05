Using barcoded templates
========================
Sequencing experiments that involve sample amplification often lose the 
information on how many different genomes/templates were actually sequenced. 
Attaching unique barcodes to the input templates before their amplification 
allows recovering the counts on the number of amplified templates. 
AccuNGS allows for long amplicons and does not assume anything about the 
length of the amplified region. Furthermore, the usage of Nextera XT to insert 
sequencing indexes and adapters at random positions results with the barcode 
in different positions for different reads (if present at all). To address this
we've written a python utility that performs an alignment of each read against 
an RT/PCR primer containing a degenerate region that allows tagging reads with 
their corresponding barcode (if present). The code is available as a python 
package here. 

