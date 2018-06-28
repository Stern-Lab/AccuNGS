from __future__ import print_function
import re, sys, argparse
from Bio import pairwise2, SeqIO, Seq

def align(seq,target):
    return pairwise2.align.localms(seq, target,10,-5,-10,-3)
    
def extract_primer_id(alignment, primer_id_start_idx, primer_id_end_idx):
    return alignment[0][alignment[3]+primer_id_start_idx:alignment[3]+primer_id_end_idx]

def primer_id_generator(fastq_handle, fasta, primer_id_start_idx, primer_id_end_idx):
    for seqRecord in SeqIO.parse(fastq_handle, "fastq"):
        toAlign=seqRecord.seq
        alignments = align(toAlign,fasta)
        if alignments[0][2] > 230:    
            yield extract_primer_id(alignments[0], primer_id_start_idx, primer_id_end_idx)
        alignments = align(toAlign.reverse_complement(),fasta)
        if alignments[0][2] > 230:        
            yield extract_primer_id(alignments[0], primer_id_start_idx, primer_id_end_idx)    

def main(args):
    #search for at least 8 N's in the provided fasta
    primer_id_template=r'(N{8,})'
    input_fasta= args.fasta
    input_fastq= args.fastq
    output_file= args.output
    record = SeqIO.read(input_fasta, "fasta")
    
    m = re.search(primer_id_template, str(record.seq))
    #reduce anchors to max of 20 bases from each side of the barcode
    recZone = record.seq[max(0,m.start(0)-20):min(len(record.seq)-1,m.end(0)+20)]
    m = re.search(primer_id_template, str(recZone))
    
    with open(input_fastq, "r") as fastq_handle:
        primer_ids = primer_id_generator(fastq_handle, recZone, m.start(0), m.end(0))
        with open(output_file, "w") as output_handle:
            for primer_id in primer_ids:
                print(primer_id, file=output_handle)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta",type=str, help="input barcode file (fasta format)")
    parser.add_argument("fastq", type=str, help="input reads file (fastq format)")
    parser.add_argument("-o","--output", type=str, help="identified barcodes file", default="identified_barcodes.txt", required=False)
    
    args = parser.parse_args(sys.argv[1:])
    main(args)