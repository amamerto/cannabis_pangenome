import pandas as pd
import sys
import warnings
warnings.filterwarnings('ignore')

def getDIR(df):
    for index, row in df.iterrows():
        if row['tstart'] < row['tstop']:
            df.loc[index, 'direction'] = '+'
            df.loc[index, 'start'] = row['tstart']
            df.loc[index, 'stop'] = row['tstop']
        else:
            df.loc[index, 'direction'] = '-'
            df.loc[index, 'start'] = row['tstop']
            df.loc[index, 'stop'] = row['tstart']
    df['start'] = df['start'].astype(int)
    df['stop'] = df['stop'].astype(int)
    return(df)

def removeDUPES(df):
    startDict = {}
    for index, row in df.iterrows():
        if row['start'] in startDict:    
            startDict[row['start']].append(index)
        else:
            startDict[row['start']] = []
            startDict[row['start']].append(index)
    for key in startDict.keys():
        if len(startDict[key]) < 2:
            pass
        else:
            keep = "none"
            for index in startDict[key]:
                if keep == "none":
                    keep = index                        
                elif df.loc[index, 'bitscore'] > df.loc[keep, 'bitscore']:
                    df.drop(keep, inplace=True)
                    keep = index
                else:
                    df.drop(index, inplace=True)

    endDict = {}
    for index, row in df.iterrows():
        if row['stop'] in endDict:    
            endDict[row['stop']].append(index)
        else:
            endDict[row['stop']] = []
            endDict[row['stop']].append(index)
    for key in endDict.keys():
        if len(endDict[key]) < 2:
            pass
        else:
            keep = "none"
            for index in endDict[key]:
                if keep == "none":
                    keep = index
                elif df.loc[index, 'bitscore'] > df.loc[keep, 'bitscore']:
                    df.drop(keep, inplace=True)
                    keep = index
                else:
                    df.drop(index, inplace=True)      
    return(df)

def cleanDF(df):
    df = removeDUPES(df)
    df = df.sort_values(by=['start'])
    df = df.reset_index(drop=True)
    
    col = ['query','target','length','direction','start','stop']
    reducedDF = pd.DataFrame(columns=col)
    cleanedDF = pd.DataFrame(columns=['query','target','ident','length','mismatch','gaps','qstart','qstop','tstart','tstop','eval','bitscore','direction','start','stop','check'])
    for gene in df['query'].unique():
        geneDF = df[df['query'] == gene]
        neg = geneDF[geneDF['direction'] == '-']
        if len(neg) > 1:
            neg = neg.iloc[::-1]
            neg = neg.reset_index(drop=True)
            neg['check'] = ""
            prev_stop = -1000000
            position = 5000
            count = -1
            check = 0

            for index,row in neg.iterrows():
                if count == -1 and row['qstart']>10:
                    pass
                elif neg.loc[index, 'qstart'] < position:
                    if count == -1: #start of table
                        count = 1
                        temp = pd.DataFrame(columns=col)
                        temp.loc[count,'query'] = row['query']
                        temp.loc[count,'target'] = row['target']
                        temp.loc[count,'stop'] = row['stop']
                        
                        position = neg.loc[index, 'qstart']
                        prev_stop = row['stop']
                        check = 1
                        neg.loc[index,'check'] = check
                        # print('START GENE')
                    elif index == len(neg)-1: #new gene end of table
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])

                        check = 1
                        neg.loc[index,'check'] = check
                        neg = neg.drop([index])
                        # print('NEW GENE - LAST ROW')
                    else: # new gene
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])

                        temp = pd.DataFrame(columns=col)
                        count = 1
                        temp.loc[count,'query'] = row['query']
                        temp.loc[count,'target'] = row['target']
                        temp.loc[count,'stop'] = row['stop']
                        position = neg.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        neg.loc[index,'check'] = check
                        # print('NEW GENE')
                elif row['start']-prev_stop > 1000000: #large gap
                    if index == len(neg)-1: #new gene end of table
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])

                        check = 1
                        neg.loc[index,'check'] = check
                        neg = neg.drop([index])
                        # print('NEW GENE - LAST ROW')
                    else: #new gene
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])

                        temp = pd.DataFrame(columns=col)
                        count += 1
                        temp.loc[count,'query'] = row['query']
                        temp.loc[count,'target'] = row['target']
                        temp.loc[count,'stop'] = row['stop']
                        position = neg.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        neg.loc[index,'check'] = check
                        # print('NEW GENE')
                else: #next exon
                    temp.loc[count,'start'] = row['start']
                    position = neg.loc[index, 'qstart']
                    prev_stop = row['stop']
                    
                    check += 1
                    neg.loc[index,'check'] = check
                    if index == len(neg)-1: #end of table
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])
                        # print('LAST ROW')
        
            cleanedDF = pd.concat([cleanedDF,neg])
            
        pos = geneDF[geneDF['direction'] == '+']
        if len(pos) > 1:
            pos = pos.reset_index(drop=True)
            pos['check'] = ""
            prev_stop = -1000000
            position = 5000
            count = -1
            check = 0

            for index,row in pos.iterrows():
                if count == -1 and row['qstart']>10:
                    pass
                elif pos.loc[index, 'qstart'] < position:
                    if count == -1: #start of table
                        count = 1
                        temp = pd.DataFrame(columns=col)
                        temp.loc[count,'query'] = row['query']
                        temp.loc[count,'target'] = row['target']
                        temp.loc[count,'start'] = row['start']
                        
                        position = pos.loc[index, 'qstart']
                        prev_stop = row['stop']
                        check = 1
                        pos.loc[index,'check'] = check
                        # print('START GENE')
                    elif index == len(pos)-1: #new gene end of table
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])
                        
                        check = 1
                        pos.loc[index,'check'] = check
                        pos = pos.drop([index])
                        # print('NEW GENE - LAST ROW')
                    else: # new gene
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])

                        temp = pd.DataFrame(columns=col)
                        count = 1
                        temp.loc[count,'query'] = row['query']
                        temp.loc[count,'target'] = row['target']
                        temp.loc[count,'start'] = row['start']
                        position = pos.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        pos.loc[index,'check'] = check
                        # print('NEW GENE')
                elif row['start']-prev_stop > 1000000: #large gap
                    if index == len(pos)-1: #new gene end of table
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])

                        check = 1
                        pos.loc[index,'check'] = check
                        pos = pos.drop([index])
                        # print('NEW GENE - LAST ROW')
                    else: #new gene
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])
                        
                        temp = pd.DataFrame(columns=col)
                        count +=1
                        temp.loc[count,'query'] = row['query']
                        temp.loc[count,'target'] = row['target']
                        temp.loc[count,'start'] = row['start']
                        position = pos.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        pos.loc[index,'check'] = check
                        # print('NEW GENE')
                else: #next exon
                    temp.loc[count,'stop'] = row['stop']
                    position = pos.loc[index, 'qstart']
                    prev_stop = row['stop']
                    
                    check += 1
                    pos.loc[index,'check'] = check
                    if index == len(pos)-1: #end of table
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        reducedDF = pd.concat([reducedDF,temp])
                        # print('LAST ROW')

            cleanedDF = pd.concat([cleanedDF,pos])
    cleanedDF = cleanedDF.sort_values(by=['start'])
    cleanedDF = cleanedDF.reset_index(drop=True)
    reducedDF = reducedDF.sort_values(by=['start'])
    reducedDF = reducedDF.reset_index(drop=True)
    return(reducedDF,cleanedDF)

def main():
    name = sys.argv[1]
    countDF = pd.DataFrame(columns=['genome','THCAS','CBDAS','CBCAS','AAE1','OAC','OLS','PT4','GPPS_ls','GPPS_ss'])
    countDF.loc[0,'genome']=name

    geneList = ['AAE1','OAC','OLS','PT4','GPPS_ls','GPPS_ss']
    for gene in geneList:
        blastfile = name+'.'+gene+'.tblastn.out'
        df = pd.read_csv(blastfile,sep='\t',names=['query','target','identity','length','mismatch','gaps','qstart','qstop','tstart','tstop','eval','bitcore'])
        chrs = list(df['target'].unique())
        count = 0
        for c in chrs:
            temp = getDIR(df[df['target']==c])
            temp = temp[(temp['identity']>75)]
            if len(temp)>0:
                Rtemp, Ctemp = cleanDF(temp)
                count = count+len(Rtemp)
            elif len(temp)==1:
                count = count+len(temp)
        countDF.loc[0,gene]=count

    synDF = pd.read_csv(name+'_filterhits.csv')
    synDF = synDF[(synDF['type']=='Full') & (synDF['length']==1635)]

    synList = ['THCAS','CBDAS','CBCAS']
    for gene in synList:
        count = len(synDF[synDF['query']==gene])
        countDF.loc[0,gene]=count

    countDF.to_csv('pathway_copynumbers.csv',index=False,header=None,mode='a')

if __name__ == '__main__':
    main()
