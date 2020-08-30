import os
import argparse
import pandas as pd
from matplotlib import pyplot as plt
from plotting import set_plots_size_params

# TODO: move to haplotypes


def plot_stretches_deep_dive(df, stretches, output_folder):
    for stretch in stretches.index:
        stretch_df = df[df.Stretch == stretch]
        stretch_dist = round(stretch_df.iloc[0, 4], 3)
        plt.figure(figsize=(20, 10))
        stretch_df.apply(lambda row: plt.plot((row['Pos1'], row['Pos2']), (row['Freq'], row['Freq']), 'C0'), axis=1)
        plt.xlabel('Position')
        plt.ylabel('Frequency')
        plt.title(f"Stretch {stretch} with meandist of {stretch_dist} made of {stretches[stretch]/2} pairs")
        plt.savefig(os.path.join(output_folder, f'stretch_{stretch}.png'))


def plot_stretches_summary(df, stretches, output_folder):
    plt.figure(figsize=(20, 10))
    for strech in stretches.index:
        strech_df = df[df.Stretch == strech]
        max_pos = max(strech_df.Pos2.max(), strech_df.Pos1.max())
        min_pos = min(strech_df.Pos1.min(), strech_df.Pos2.min())
        meandist = round(strech_df.iloc[0, 4], 3)
        plt.plot((min_pos, max_pos), (meandist, meandist), label=f"Stretch #{strech}")
    plt.legend()
    plt.title("Stretches")
    plt.xlabel('Position')
    plt.ylabel('Mean Distance')
    plt.savefig(os.path.join(output_folder, f'haplotype_overview.png'))


def graph_haplotypes(input_file, number_of_stretches, output_dir):
    set_plots_size_params(20)
    df = pd.read_csv(input_file, sep="\t")
    biggest_stretches = df.Stretch.value_counts()[:number_of_stretches]
    plot_stretches_deep_dive(df=df, stretches=biggest_stretches, output_folder=output_dir)
    plot_stretches_summary(df=df, stretches=biggest_stretches, output_folder=output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", help="Where all the files will go", required=True)
    parser.add_argument("-i", "--input_file", help="stretches.csv file created by co_occurs_to_stretches.py",
                        required=True)
    parser.add_argument("-n", "--number_of_stretches", help="plot this many stretches, default is 5", default=5,
                        type=int)
    args = parser.parse_args()
    graph_haplotypes(input_file=args.input_file, number_of_stretches=args.number_of_stretches,
                     output_dir=args.output_dir)
