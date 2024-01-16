import pandas as pd
import numpy as np
from plotnine import *
from argparse import ArgumentParser
from pandas.api.types import CategoricalDtype
import seaborn as sns

def parse_arguments():
    parser = ArgumentParser(description='pandots')
    parser.add_argument(
        '--paf',
        required=True,
        help='paf alignment file'
    )
    parser.add_argument(
        '--ref',
        default='Reference',
        help='reference chromosome name'
    )
    parser.add_argument(
        '--out',
        default='pandots',
        help='outfile prefix name'
    )
    parser.add_argument(
        '--pad',
        default=10000000,
        help='pad value between genome lines'
    )
    parser.add_argument(
        '--key',
        type=bool,
        choices=[True, False],
        default=False,
        help='show genome key or not'
    )
    parser.add_argument(
        '--flip',
        type=bool,
        choices=[True, False],
        default=False,
        help='make paf column 1 ref or query'
    )
    parser.add_argument(
        '--invert',
        type=bool,
        choices=[True, False],
        default=False,
        help='invert reference in relation to pangenome'
    )
    parser.add_argument(
        '--reorder',
        type=bool,
        choices=[True, False],
        default=False,
        help='adjust order of genomes plotted'
    )
    parser.add_argument(
        '--recolor',
        type=bool,
        choices=[True, False],
        default=False,
        help='adjust color of genomes plotted'
    )
    return parser.parse_args()

def invert_paf(df, orient):
    for genome in df['queryID'].unique():
        if orient.get(genome, False):
            subset = df[df['queryID'] == genome]

            newstart = [l - stop for l, stop in zip(subset['queryLen'], subset['queryEnd'])]
            newend = [l - start for l, start in zip(subset['queryLen'], subset['queryStart'])]
            newstrand = ["-" if s == "+" else "+" for s in subset['strand']]

            df.loc[subset.index, 'queryStart'] = newstart
            df.loc[subset.index, 'queryEnd'] = newend
            df.loc[subset.index, 'strand'] = newstrand
    return df

def parital_invert_paf(df, orient):
    for genome in df['queryID'].unique():
        if orient.get(genome, False):
            subset = df[df['queryID'] == genome]

            newstart = [start if s == '+' else stop for s, start, stop in zip(subset['strand'], subset['queryStart'], subset['queryEnd'])]
            newend = [stop if s == '+' else start for s, start, stop in zip(subset['strand'], subset['queryStart'], subset['queryEnd'])]

            df.loc[subset.index, 'queryStart'] = newstart
            df.loc[subset.index, 'queryEnd'] = newend
    return df

def load_paf(fpath):
    df = pd.read_csv(fpath, sep='\t', header=None, names=range(40))
    df = df.iloc[:, :12]
    df.columns = [
        "refID","refLen","refStart","refEnd","strand",
        "queryID","queryLen","queryStart","queryEnd",
        "numResidueMatches","lenAln","mapQ"
    ]
    return df

def load_orients(fpath):
    df = pd.read_csv(fpath,sep='\t')
    orients = df.set_index('Chromosome').to_dict()
    orients = orients['Flip']
    return orients

def stack_query(df,pad):
    add = 0
    for genome in df['queryID'].unique():
        df['queryStart'] = df.apply(lambda x: x['queryStart']+add if x['queryID']==genome else x['queryStart'], axis=1)
        df['queryEnd'] = df.apply(lambda x: x['queryEnd']+add if x['queryID']==genome else x['queryEnd'], axis=1)
        add = add + pad
    return df

def draw_ggplot(pairs, xlab, ylab, tlab, klab, key, recolor):
    df = pd.DataFrame([(x1, y1, x2, y2, genome) for (x1, y1), (x2, y2), genome in pairs],
                      columns=['xstart', 'ystart', 'xstop', 'ystop', 'genome'])
    if not recolor:
        gp = ggplot(df, aes(x='xstart', y='ystart', xend='xstop', yend='ystop', color='genome'))
    else:
        df['color'] = df['genome'].astype(CategoricalDtype(categories=df['genome'].unique(), ordered=True))
        gp = ggplot(df, aes(x='xstart', y='ystart', xend='xstop', yend='ystop', color='color'))
        viridis_palette = sns.color_palette("viridis", as_cmap=True)
        new_viridis_palette = viridis_palette(np.linspace(0, 1, 256))[1:]
        color_palette = new_viridis_palette
        gp += scale_color_manual(values=sns.color_palette(color_palette, n_colors=len(df['genome'].unique())).as_hex())
    gp = gp + geom_segment() + \
        labs(x=xlab, y=ylab, title=tlab, color=klab) + \
        scale_x_continuous(labels=lambda breaks: [f'{int(b/1e6)}Mb' for b in breaks]) + \
        scale_y_continuous(labels=lambda breaks: [f'{int(b/1e6)}Mb' for b in breaks]) + \
        theme_minimal() + \
        theme(axis_text_y=element_blank(), figure_size=(12, 8))
    if not key:
        gp += theme(legend_position='none', axis_text_y=element_blank())
    return gp

def reorder_genomes(df):
    df['queryAlnLen'] = df['queryEnd'] - df['queryStart']
    df['queryContigLen'] = df.groupby('queryID')['queryAlnLen'].transform('sum')
    neg_vals = df[df['queryAlnLen'] < 0].groupby('queryID')['queryAlnLen'].sum().reset_index(name='neg_vals')
    df = pd.merge(df, neg_vals, on='queryID', how='left')
    df = df.sort_values(by=['queryContigLen'])
    return df

def main():
    args = parse_arguments()

    paf = pd.read_csv(args.paf, sep='\t', header=None, names=range(40))
    paf = paf.iloc[:, :12]
    if args.flip:
        paf.columns = [
            "queryID","queryLen","queryStart","queryEnd","strand",
            "refID","refLen","refStart","refEnd",
            "numResidueMatches","lenAln","mapQ"
        ]
    else:
        paf.columns = [
            "refID","refLen","refStart","refEnd","strand",
            "queryID","queryLen","queryStart","queryEnd",
            "numResidueMatches","lenAln","mapQ"
        ]
    orients = load_orients("csat_orientations.tsv")
    if args.invert:
        orients = {k :not b for k,b in orients.items()}

    paf = invert_paf(paf,orients)
    paf = parital_invert_paf(paf,orients)
    if args.reorder:
        paf = reorder_genomes(paf)
    paf = stack_query(paf,args.pad)

    paf['queryName'] = paf['queryID'].str.split('.').str[0]
    pairs = [
        ((a, b), (c, d), e)
        for a,b,c,d,e in zip(paf.refStart, paf.queryStart, paf.refEnd, paf.queryEnd, paf.queryName)
    ]
        
    plot = draw_ggplot(pairs,args.ref,'','','Genome',args.key,args.recolor)

    ggsave(plot,args.out+'.pandots.svg',format="svg")


if __name__ == '__main__':
    main()