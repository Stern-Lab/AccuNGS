import sys, gzip, argparse
from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord

def merger_generator(forward_handle,reverse_handle, rep_length):
    for a, b in zip (SeqIO.parse(forward_handle, "fastq"),SeqIO.parse(reverse_handle, "fastq")):
        if (a.id.split(" ")[0] !=  b.id.split(" ")[0]):
            print("Problem:" + str(i))
        new_seq_id=a.id.split(" ")[0]
        new_seq_str = str(a.seq) + ("N"*rep_length) + str(b.seq)
        a_quals=a.letter_annotations["phred_quality"]
        b_quals=b.letter_annotations["phred_quality"]
        new_seq_qual= a_quals+[1.0 for a in range(rep_length)]+b_quals

        new_seq=SeqRecord(Seq.Seq(new_seq_str),id=new_seq_id,description="",letter_annotations={"phred_quality":new_seq_qual})
        yield new_seq

def main(args):
    file1 = args.forward
    file2 = args.reverse
    output_file = args.output_file
    rep_length = args.rep_length

    if file1.endswith(".gz"):
        concat_fastq_gz(file1, file2, output_file, rep_length)
    else:
        concat_fastq(file1, file2, output_file, rep_length)

def concat_fastq(file1, file2, output_file, rep_length):
    with open(file1, "r") as forward_handle:
        with open(file2, "r") as reverse_handle :
            with open (output_file,"w") as merged_handle:
                SeqIO.write(merger_generator(forward_handle,reverse_handle, rep_length), merged_handle, "fastq")

def concat_fastq_gz(file1, file2, output_file, rep_length):
    with gzip.open(file1, "rt") as forward_handle:
        with gzip.open(file2, "rt") as reverse_handle :
            with gzip.open (output_file,"wt") as merged_handle:
                SeqIO.write(merger_generator(forward_handle,reverse_handle, rep_length), merged_handle, "fastq")
                
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("forward",type=str, help="Forward FASTQ file")
    parser.add_argument("reverse",type=str, help="Reverse FASTQ file")
    parser.add_argument("output_file", type=str, help="output: merged FASTQ file")
    parser.add_argument("-n", "--rep_length", help="amount of N bases to repeat", type=int, default=60, required=False)
    
    args = parser.parse_args(sys.argv[1:])
    main(args)

