import sys, argparse
import numpy as np
import pandas as pd


class comutations_bubble(object):
    def __init__(self, a, b, freq, pval):
        self.nodes=set([a,b])
        self.distances={'_'.join([str(a),str(b)]):freq, '_'.join([str(b),str(a)]):freq}
        self.pvalues={'_'.join([str(a),str(b)]):pval, '_'.join([str(b),str(a)]):pval}
        self.meandist=np.mean(list(self.distances.values()))

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        return ':'.join([str(a) for a in [len(self),min(self.nodes),max(self.nodes)]])

    def can_merge(self, other, distance=10):
        return len(self.nodes.intersection(other.nodes)) > 0 and ((other.meandist*distance >= self.meandist) and (other.meandist/distance <= self.meandist))
    
    def union(self, other):
        self.nodes=self.nodes.union(other.nodes)
        self.distances.update(other.distances)
        self.pvalues.update(other.pvalues)
        self.meandist=np.mean(list(self.distances.values()))


def load_file(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = ["start", 'start_mutation', "end", "end_mutation", "pval", "variant_freq"]
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    df = df[df["pval"]<1]
    return df


def obtain_comutations(comutations, max_pval=10**-9, distance=10):
    dfs=[]
    sig_positions = comutations[(comutations["pval"] < max_pval)]
    lines=sig_positions.itertuples(index=False)
    nodes=[]
    for line in lines:
        a = line[0]
        b = line[2]
        pval = line[4]
        dist = line[5]
        node = comutations_bubble(a, b, dist, pval)
        nodes.append(node)

    results=[]
    nodes=sorted(nodes, key=lambda item: -item.meandist)
    while len(nodes) > 0:
        bubble=nodes.pop(0)
        merged=False
        for item in nodes:
            if item.can_merge(bubble, distance):
                item.union(bubble)
                merged=True
                break
        if not merged:
            results.append(bubble)
                
    if results:
        #Collate lines
        for i, cluster in zip(range(len(results)), sorted(results, key=lambda item: -item.meandist)):
            distances=map(lambda x: (int(x.split("_")[0]), int(x.split("_")[1]), cluster.distances[x]), cluster.distances.keys())           
            data = pd.DataFrame(distances, columns=["Pos1", "Pos2", "Freq"])
            data["Stretch"]=i

            data["meandist"]=cluster.meandist
            dfs.append(data)
    
    if len(dfs) == 0:
        return pd.DataFrame({"Pos1":[], "Stretch":[], "meandist":[]})    
    return pd.concat(dfs) 


def calculate_stretches(linked_pairs, max_pval, distance, output):
    if not max_pval:
        max_pval = 10 ** -9
    if distance is None:
        distance = 10
    comutations = load_file(linked_pairs)
    res_table = obtain_comutations(comutations, max_pval, distance)
    kws={"index":False, "sep":"\t"}
    if not output:
        res_table.to_csv(linked_pairs + ".stretches.out", **kws)
    else:
        res_table.to_csv(output, **kws)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "linked_pairs", type=str, help="all linked pairs to find haplotypes"
    )
    parser.add_argument(
        "-p", "--max_pval",
        type=float,
        help="Maximal p-value of a given pair to consider for merging with another pair",
        default=10 ** -9,
    )
    parser.add_argument(
        "-d", "--distance",
        type=float,
        help="The position to consider pairs",
        default=10,
    )
    parser.add_argument(
        "-o", "--out",
        type=str,
        help="The output file",
        default=None,
    )

    args = parser.parse_args(sys.argv[1:])
    calculate_stretches(linked_pairs=args.linked_pairs, distance=args.distance, max_pval=args.max_pval, output=args.out)
