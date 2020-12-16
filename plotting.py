import os

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def set_plots_size_params(size):
    # Adapted from https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
    bigger = size * 1.2
    slightly_bigger = size * 1.1
    plt.rc('font', size=size)                        # controls default text sizes
    plt.rc('axes', titlesize=bigger)                 # fontsize of the axes title
    plt.rc('axes', labelsize=slightly_bigger)        # fontsize of the x and y labels
    plt.rc('xtick', labelsize=size)                  # fontsize of the tick labels
    plt.rc('ytick', labelsize=size)                  # fontsize of the tick labels
    plt.rc('legend', fontsize=size)                  # legend fontsize


def graph_read_counter(read_counter_file, ax):
    data = pd.read_csv(read_counter_file, sep="\t")
    data.number_of_alignments.hist(bins=10, ax=ax)
    ax.grid(False)
    ax.set_ylabel("Density")
    ax.set_xlabel("Times mapped")
    ax.title.set_text('Number of Reads')
    return ax


def graph_blast_length_distribution(blast_file, ax):
    data = pd.read_csv(blast_file, sep="\t", usecols=[7])
    data.hist(bins=30, ax=ax)
    ax.grid(False)
    ax.set_ylabel("Density")
    ax.set_xlabel("Alignment Length")
    ax.title.set_text('Length distribution')
    return ax


def graph_coverage(freqs, ax):
    data = freqs[freqs.base_rank == 0]
    ax.plot(data.ref_pos, data.coverage)
    ax.set_ylabel("Coverage")
    ax.set_xlabel("Reference Position")
    ax.title.set_text('Coverage')
    return ax


def graph_mutation_freqs_by_mutation(mutation_data, ax):
    mutation_data["mutation"] = mutation_data["ref_base"] + ">" + mutation_data["read_base"]
    sns.boxplot(x='mutation', y='frequency', data=mutation_data, ax=ax)
    ax.set_yscale("log")
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Mutation")
    ax.title.set_text('Mutations by type')
    return ax


def graph_mutation_freqs_by_position(mutation_data, ax):
    sns.scatterplot(x="ref_pos", y="frequency", data=mutation_data, ax=ax)
    ax.set_ylim([mutation_data['frequency'].min() * 0.9, 1])  # otherwise yaxis is missing lower values.. weird.
    ax.set_yscale("log")
    ax.set_xlabel("Ref position")
    ax.set_ylabel("Frequency")
    ax.title.set_text(f"Mutations frequencies by position (sum of mutations: {len(mutation_data)})")
    return ax


def plot_stretches_summary(df, stretches, ax):
    plt.figure(figsize=(20, 10))
    for strech in stretches.index:
        strech_df = df[df.Stretch == strech]
        max_pos = max(strech_df.Pos2.max(), strech_df.Pos1.max())
        min_pos = min(strech_df.Pos1.min(), strech_df.Pos2.min())
        meandist = round(strech_df.iloc[0, 4], 3)
        ax.plot((min_pos, max_pos), (meandist, meandist), label=f"Stretch #{strech}")
    if not df.empty:
        ax.legend() # TODO : add to log if df is empty and not preformed
    ax.title.set_text("Haplotypes")
    ax.set_xlabel('Position')
    ax.set_ylabel('Mean Frequency')
    return ax


def graph_haplotype_overview(stretches_file, ax, stretches_to_plot):
    if not os.path.isfile(stretches_file):
        ax.text(0.5, 0.5, "Waiting for data...", fontsize=18, ha='center')
        ax.title.set_text("Haplotypes")
        return ax
    df = pd.read_csv(stretches_file, sep="\t")
    biggest_stretches = df.Stretch.value_counts()[:stretches_to_plot]
    return plot_stretches_summary(df=df, stretches=biggest_stretches, ax=ax)


def graph_summary(freqs_file, blast_file, read_counter_file, stretches_file, output_file, stretches_to_plot,
                  min_coverage):
    # TODO: graph multimapped ignored bases, min_coverage can be set as defult
    if not min_coverage:
        min_coverage = 10
    fig, axes = plt.subplots(figsize=(60, 20), ncols=3, nrows=2)
    plt.suptitle("Pipeline Statistics", fontsize=18)
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    freqs = pd.read_csv(freqs_file, sep="\t")
    mutation_data = freqs[(freqs['base_rank'] != 0) & (freqs['coverage'] > min_coverage) &
                          (freqs['base_count'] > 0)].copy()  # TODO: in the old pipeline theres also a probability thing here
    axes[0][0] = graph_blast_length_distribution(blast_file, axes[0][0])
    axes[0][1] = graph_read_counter(read_counter_file, axes[0][1])
    axes[0][2] = graph_coverage(freqs, axes[0][2])
    axes[1][0] = graph_mutation_freqs_by_mutation(mutation_data, axes[1][0])
    axes[1][1] = graph_mutation_freqs_by_position(mutation_data, axes[1][1])
    axes[1][2] = graph_haplotype_overview(stretches_file, axes[1][2], stretches_to_plot=stretches_to_plot)
    fig.savefig(output_file, bbox_inches="tight", pad_inches=0.3)
