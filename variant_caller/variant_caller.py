import sys, argparse
import pandas as pd
import numpy as np
from scipy import stats

def learn_gammas(data):
    """
    Fits gamma distributions on control to several mutation types and possibly their broader dinucleotide context
    """
    gamma_distributions=[]
    grouped = data[(data["Rank"] != 0)].groupby(["mutation_type"], as_index=False)
    for name, group in grouped:
        shape, loc, scale= fit_gamma_iterative(name, group)
        gamma_distributions.append([name, "all", shape, scale, stats.gamma.ppf(0.95, shape, loc, scale), group["adjusted_freq"].max(), group["adjusted_freq"].mean(), group["adjusted_freq"].median()])

    #validate gamma on transition mutations
    grouped = data[(data["Rank"] != 0) & (data["mutation_type"] == "CT")].groupby(["next"], as_index=False)
    for name, group in grouped:
        shape, loc, scale= fit_gamma_iterative(name, group)
        gamma_distributions.append(["CT","X" + name, shape, scale, stats.gamma.ppf(0.95, shape, loc, scale), group["adjusted_freq"].max(), group["adjusted_freq"].mean(), group["adjusted_freq"].median()])
        
    grouped = data[(data["Rank"] != 0) & (data["mutation_type"] == "GA")].groupby(["prev"], as_index=False)
    for name, group in grouped:
        shape, loc, scale= fit_gamma_iterative(name, group)
        gamma_distributions.append(["GA",name+"X", shape, scale, stats.gamma.ppf(0.95, shape, loc, scale), group["adjusted_freq"].max(), group["adjusted_freq"].mean(), group["adjusted_freq"].median()])       
    
    return gamma_distributions

def fit_gamma_iterative(name, group):
    """
    Iteratively fits a gamma distribution to the provided group.
    Outliers are identified and removed after fitting by the three-sigma-rule, and then refit
    """
    converged=False
    tofit=group
    while not converged:
        param = stats.gamma.fit(tofit["adjusted_freq"]  , floc=0)
        std=stats.gamma.std(param[0],param[1],param[2])
        mean=stats.gamma.mean(param[0],param[1],param[2])
        max_obs=tofit["adjusted_freq"].max()
        if (mean+3*std) < max_obs: #three-sigma-rule for identiftying outliers
            tofit = tofit[tofit["adjusted_freq"] != max_obs] #exclude extreme outlier and reiterate
        else:
            return param[0],param[1],param[2]

def call_by_gamma(freq, mutation_type, previous, forward, gammas):
    """
    :return: the probability for a given frequency to arise from the corresponding gamma distribution
    """
    matching_lines = [x for x in gammas if x[0] == mutation_type]
    if len(matching_lines) > 1:
        if mutation_type == "GA":
            matching_lines = [x for x in matching_lines if x[1] == str(previous+"X")]
        if mutation_type == "CT":
            matching_lines = [x for x in matching_lines if x[1] == str("X"+forward)]
    if len(matching_lines) == 0:
        return 1
    line = matching_lines[0]
    p_val = stats.gamma.cdf(freq,float(line[2]),0,float(line[3]))
    return 1-p_val

def output_to_file(data, output_file):
    """
    Outputs a file similar in structure to the output file of V-Phaser 2.0 by Broad (Yang 2013, BMC Genomics).
    Indels are currently not supported     
    """

    dataFiltered=data[(data["Rank"]!=0)]
    
    dataToOutput=pd.DataFrame()
    dataToOutput["Ref_Pos"]=dataFiltered["Pos"]
    dataToOutput["Var"]=dataFiltered["Base"]
    dataToOutput["Cons"]=dataFiltered["Ref"]   
    dataToOutput["pval"]=dataFiltered["to_call"]
    dataToOutput["Type"]="snp"
    dataToOutput["Var_perc"]=dataFiltered["Freq"]
    dataToOutput["SNP_Profile"]=dataFiltered["Ref"] + ":" + dataFiltered["major_read_count"].astype(np.int32).astype(str)  + " " + dataFiltered["Base"] + ":" + dataFiltered["counts_for_position_x"].astype(np.int32).astype(str) 
    dataToOutput["pval"] = np.where(dataToOutput["Var_perc"] == 0, "1.0", dataToOutput["pval"])
    dataToOutput["pval"]=dataToOutput["pval"].apply(pd.to_numeric)

    print ("writing output to " + output_file)
    dataToOutput.to_csv(output_file, index=False)

def main(args):
    """
    Identifies the probabilities of variant frequencies to come from their corresponding background errors distributions
    
    :param Args[0] - sample file
    :param Args[1] - homogeneous control file
    :param Args[2] - minimum coverage to consider positions when fitting distributions and identifying variants 
    
    :return: "output.var.csv"
    """
    sample_file = args.sample
    control_file = args.control
    min_coverage = args.coverage

    label_control = "Control"
    label_sample = "Sample"

    print ("loading " + sample_file + " as sample")
    data_mutations = pd.read_table(sample_file)
    data_mutations["source"] = label_sample
    print ("loading " + control_file + " as homogeneous control")
    data_control = pd.read_table(control_file)
    data_control["source"] = label_control
    
    data = pd.concat([data_control, data_mutations])
    #generate consensus of both sample and control
    consensus_data = data[['Pos','Base','source',"Read_count","Freq"]][data['Rank'] == 0]
    consensus_data = consensus_data[consensus_data['Pos'] == np.round(consensus_data['Pos'])] #removes insertions
    consensus_data['next'] = consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data['prev'] = consensus_data.shift(1)['Base'] + consensus_data['Base'] 
    consensus_data['Pos'] = consensus_data[['Pos']].apply(pd.to_numeric)
    consensus_data["context"] = consensus_data.shift(1)['Base'] + consensus_data['Base'] + consensus_data.shift(-1)['Base']
    consensus_data["major_read_count"] = np.round(consensus_data['Read_count']*consensus_data['Freq'])
    consensus_data = consensus_data.rename(columns={"Base":"consensus"})
    data = pd.merge(data, consensus_data[['Pos','prev','next','context','consensus',"source","major_read_count"]], on=["Pos","source"])    
    
    data['counts_for_position'] = np.round(data['Read_count']*data['Freq'])
    data = data[data['Pos'] == np.round(data['Pos'])] #remove insertions
    data['Pos'] = data[['Pos']].apply(pd.to_numeric)
    data = data[(data['Read_count'] > min_coverage)]
    data['mutation_type'] = data['prev'].str[1] + data['Base']
    data['mutation_class'] = np.where(data["Rank"] == 0,"self",np.where(data["mutation_type"].str.contains("-"),"deletion",
                                                                        np.where(data['mutation_type'].isin(['GA','AG','CT','TC']), 'transition', 'transversion')))
    data = data[data["mutation_class"] != "deletion"]
    grouped = data.groupby(["source","Pos"], as_index=False)["counts_for_position"].agg("sum")
    data = pd.merge(data,grouped, on=["source","Pos"])
    data["adjusted_freq"] = (np.where(1>data["counts_for_position_x"],1,data["counts_for_position_x"]) / data["counts_for_position_y"])   
       
    gamma_distributions=[]
    gamma_distributions.extend(learn_gammas(data[data["source"] == label_control]))
    
    to_call = data[data["source"] == label_sample].copy()
    to_call.loc[:,"to_call"] = to_call.apply(lambda row: call_by_gamma(row["adjusted_freq"],row["mutation_type"], row["prev"],row["next"],gamma_distributions), axis=1)

    output_to_file(to_call[(to_call["source"]==label_sample)], args.output)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sample",type=str, help="input sample .freqs file")
    parser.add_argument("control", type=str, help="input control .freqs file")
    parser.add_argument("-c", "--coverage", type=int, help="minimum position coverage to fit and call variants", required=False, default=100000)
    parser.add_argument("-o", "--output", type=str, help="output variant file", required=False, default="output.var.csv")
    
    args = parser.parse_args(sys.argv[1:])
    main(args)
