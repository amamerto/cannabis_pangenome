import pandas as pd
import sys
import warnings
warnings.filterwarnings('ignore')

def main():
    tbl = sys.argv[1]
    name = tbl.split('_')[0]

    df = pd.read_csv(tbl)
    df = df[(df['type']=='Full') & (df['length']==1635)]
    df['count'] = df.index
    df['count'] = df['count'].astype(str)
    df['name'] = df['target']+'_'+df['query']+'_'+df['count']
    df['score']=0

    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)

    finaldf = df[['target','start','end','name','score','direction']]
    finaldf.to_csv(name+'_fullsyn.bed',index=False,header=False,sep='\t')

if __name__ == '__main__':
    main()