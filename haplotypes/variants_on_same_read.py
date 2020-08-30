import itertools
import sys
import argparse

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

def main(args):

    input_x=args.position

    freqs=pd.read_csv(args.freqs_file, sep="\t")
    freqs = freqs[freqs['Pos'] == np.round(freqs['Pos'])] #remove insertions
    if (input_x < freqs["Pos"].min()) or (input_x > freqs["Pos"].max()):
        sys.exit()

    all_mappings=pd.read_csv(args.blast_output, names=["read_id","start","end"], sep="\t")
    all_mutations=pd.read_csv(args.mutations_all, names=["pos","read_id","mutant","read_positions"], sep="\t")
  
    cons = freqs[(freqs["Rank"] == 0)
                 & (freqs["Base"] != "-")]
    cons.insert(0, "pos", pd.to_numeric(cons.loc[:,"Pos"]))

    all_mutations = pd.merge(all_mutations, cons[["pos","Ref"]], on="pos")
    #Identify co-occurring variants up to max overlap length - 250 bases
    variants_combinations=range(input_x+1,input_x+250)
    
    for y in variants_combinations:
        x=input_x
        maps_for_two_pos = all_mappings[(all_mappings["start"] <= x) & (all_mappings["end"] >= y)]
        grouped = maps_for_two_pos.groupby('read_id')
        grouped = grouped.filter(lambda x: len(x) == 2)
        # for every read_id which appears exactly twice (why?) merge with this specific x position
        merged = pd.merge(pd.DataFrame({"read_id":grouped["read_id"].unique()}), all_mutations[all_mutations["pos"]==x][["pos","read_id"]], on="read_id", how="left")
        # merge every mutation in position x with mutations on the same read in position y
        merged = pd.merge(merged, all_mutations[all_mutations["pos"]==y][["pos","read_id"]], on="read_id", how="left")
        x_label = "pos_" + str(x)
        y_label = "pos_" + str(y)
        merged[x_label] = np.where(merged["pos_x"] == x, 1, 0)
        merged[y_label] = np.where(merged["pos_y"] == y, 1, 0)
        ct = pd.crosstab(merged[x_label], merged[y_label])
        if ct.shape == (2,2):  # what are the odds that the mutations in x and y were observed together by chance
            fisher_test = fisher_exact(ct, alternative='greater')
            print('\t'.join([str(x) for x in [x, y, fisher_test[0], fisher_test[1], ct[1][1]*1.0/(ct[0][0]+ct[0][1]+ct[1][0]+ct[1][1])]]))
        else:
            print('\t'.join([str(x) for x in [x, y, 0.0, 1.0, 0.0]]))
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_output",type=str, help="all BLAST output for this sample")
    parser.add_argument("mutations_all", type=str, help="mutations_all.txt file (filtered from text)")
    parser.add_argument("position", type=int, help="The position to consider pairs")
    parser.add_argument("freqs_file", type=str, help="freqs file")
    
    args = parser.parse_args(sys.argv[1:])
    main(args)
