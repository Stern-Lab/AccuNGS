import sys
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
    
def load_file(path, label):
    df = pd.read_csv(path, sep="\t", names=["start","end","fisher_stat","pval","variant_freq"])
    df["sample"] = label
    df = df[df["pval"]<1]
    return df

def obtain_comutations(comutations, max_pval=10**-9, distance=10):
    dfs=[]
    sig_positions = comutations[(comutations["pval"] < max_pval)]
    
    lines=sig_positions.itertuples(index=False)
    nodes=[]
    for line in lines:
        a = line[0]
        b = line[1]
        pval = line[3]
        dist = line[4]
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

            data["Sample"]=comutations.iloc[0]['sample']
            data["meandist"]=cluster.meandist
            dfs.append(data)
    
    if len(dfs) == 0:
        return pd.DataFrame({"Pos1":[], "Stretch":[], "meandist":[]})    
    return pd.concat(dfs) 

