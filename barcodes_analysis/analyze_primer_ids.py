import re, itertools, csv, numpy, os.path, sys
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
        print ("primers file found!")
        with open(barcodes_file,"r") as handle:
            return handle.read().splitlines()  
    raise Exception("No barcodes file") 

'''
args[0] = input file with primer_ids (each primer in each line)
args[1] = output path for histogram of primer-ID abundance
'''
def main(args):
    matches = [x for x in obtain_matches(args[0]) if x.count("-") < 3]
    
    with_counts = Counter(matches) 
    
    print (with_counts)
    print ("total distinct primer IDs",len(with_counts.keys()))
    print ("total primer IDs", len(matches))
    if len (args) > 2:
        print ("total primer IDs with abundance above parameter",len(list(filterbyvalue(with_counts.values(), args[2]))))
    
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
        plt.savefig(args[1] + "histogram." + save_format , dpi=680, format=save_format )
    plt.close()

if __name__ == "__main__":
    main(sys.argv[1:])