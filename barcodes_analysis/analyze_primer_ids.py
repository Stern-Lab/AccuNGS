import re, itertools, csv, numpy, os.path, sys, argparse
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns

def filterbyvalue(seq, value):
   for el in seq:
       if el>=value: yield el

def hamming_distance(s1, s2):
    """
    :return: the Hamming distance between the two sequences of equal-length
    """
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def obtain_matches(barcodes_file):

    if os.path.isfile(barcodes_file):
        with open(barcodes_file,"r") as handle:
            return handle.read().splitlines()  
    raise Exception("No barcodes file") 

def main(args):
    matches = [x for x in obtain_matches(args.barcodes_file) if x.count("-") < 3]
    
    with_counts = Counter(matches) 
    
    print ("total primer IDs in input file: " + str(len(matches)))
    print ("total distinct primer IDs: " + str(len(with_counts.keys())))
    print ("total distinct primer IDs with abundance at or above parameter: " + str(len(list(filterbyvalue(with_counts.values(), args.cutoff)))))
    
    hist = numpy.histogram(with_counts.values(), max(with_counts.values()))
   
    xs = hist[1].tolist()
    del xs[-2]

    plt.scatter(x=xs, y=hist[0])
    plt.title("Primer IDs distribution")
    plt.xlabel("Number of raw reads per unique Primer-ID")
    plt.ylabel("Number of distinct Primer-IDs")
    plt.yscale("symlog",linthreshy=1)
    plt.xlim(1,max(with_counts.values())+10)
    plt.ylim(1,len(with_counts))
    
    for save_format in ["png"]:
        plt.tight_layout()
        plt.savefig(args.output_dir + "/histogram." + save_format , dpi=680, format=save_format )
    plt.close()

    print ("output written to " + args.output_dir + "/histogram." + save_format )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("barcodes_file",type=str, help="text file with identified barcodes")
    parser.add_argument("output_dir", type=str, help="folder for results output")
    parser.add_argument("-c","--cutoff", type=int, help="barcodes abundance cutoff", default=1, required=False)
    
    args = parser.parse_args(sys.argv[1:])
    main(args)