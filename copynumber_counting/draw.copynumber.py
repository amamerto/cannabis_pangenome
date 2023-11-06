import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def main():
    counts = sys.argv[1]

    df = pd.read_csv(counts,
                    header=0,
                    names=['genome','THCAS','CBDAS','CBCAS','AAE1','OAC','OLS','PT4']
                    )

    sns.set(rc={'figure.figsize':(15,10)})
    sns.set(font_scale=1.75)
    ax = sns.boxplot(data=df)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Counts")
    ax.set_title("Cannabinoid Pathway Counts Per Genome")
    plt.savefig("pathway_copynumber.boxplot.svg", format="svg")

if __name__ == '__main__':
    main()
