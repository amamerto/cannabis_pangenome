import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def main():
    counts = sys.argv[1]

    temp = pd.read_csv(counts,
                    header=0,
                    names=['genome','THCAS','CBDAS','CBCAS','AAE1','OAC','OLS','PT4','GPPS_ls','GPPS_ss']
                    )
    df = temp[['genome','AAE1','OLS','OAC','GPPS_ls','GPPS_ss','PT4','CBDAS','THCAS','CBCAS']]
    
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set(rc={'figure.figsize':(11,11)})
    sns.set_theme(style="ticks", palette="bright", font_scale=2, rc=custom_params)
    ax = sns.boxplot(data=df,orient="h",linewidth=3,width=0.35)
    
    ax.set_xlabel("Copy Number")
    # ax.set_ylabel("Counts")
    # ax.set_title("Cannabinoid Pathway Counts Per Genome")
    ax.yaxis.set_ticks([])  # Hide ticks
    ax.yaxis.set_ticklabels([])  # Hide tick labels

    colors = sns.color_palette("bright", n_colors=len(df.columns) - 1)  # Adjust the palette if needed
    labels = df.columns[1:]  # Exclude the 'genome' column
    legend_elements = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=30)
                    for color, label in zip(colors, labels)]
    legend = ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(2, 0, 0.35, 1), ncol=4, frameon=True, framealpha=1)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_linewidth(1.5)  # Adjust the line width
    legend.get_frame().set_edgecolor('black')  # Adjust the edge color

    plt.savefig("pathway_copynumber.boxplot.svg", format="svg")

if __name__ == '__main__':
    main()
