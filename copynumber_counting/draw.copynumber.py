import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def main():
    counts = sys.argv[1]

    temp = pd.read_csv(counts,
                    header=0,
                    names=['genome','THCAS','CBDAS','CBCAS','AAE1','OAC','OLS','PT4']
                    )
    df = temp[['genome','AAE1','OLS','OAC','PT4','CBDAS','THCAS','CBCAS']]
    
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set(rc={'figure.figsize':(10,10)})
    sns.set_theme(style="ticks", palette="tab10", font_scale=2, rc=custom_params)
    ax = sns.boxplot(data=df,orient="h",linewidth=7,width=0.35)
    ax.set_xlabel("Copy Number")
    # ax.set_ylabel("Counts")
    # ax.set_title("Cannabinoid Pathway Counts Per Genome")
    plt.savefig("pathway_copynumber.boxplot.svg", format="svg")

if __name__ == '__main__':
    main()
