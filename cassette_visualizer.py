import pandas as pd
import svgwrite as svg
from argparse import ArgumentParser
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    parser = ArgumentParser(description='cassette visualizer')
    parser.add_argument(
          '--type',
          choices=['out', 'bed', 'gff'],
          required=True,
          help='extension of input file'
    )
    parser.add_argument(
        '--file',
        required=True,
        help='input file'
    )
    parser.add_argument(
        '--scale',
        default=250,
        type=int,
        help='scale factor to size figure. ~250 for short regions. ~20000 for chromosome. default:250'
    )
    parser.add_argument(
        '--name',
        default='sample',
        type=str,
        help='prefix label to name outfiles. default:sample'
    )
    parser.add_argument(
        '--shape',
        choices=['rect','tri','flag','trap','penta','pentaI'],
        default='rect',
        type=str,
        help='shape of gene on cassette. default:rect'
    )
    parser.add_argument(
        '--pad',
        default=5,
        type=int,
        help='padding for gene shape width. default:5'

    )
    return parser.parse_args()

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

    col = ['chr','start','stop','gene','length','direction']
    exonDF = pd.DataFrame(columns=col)
    intronDF = pd.DataFrame(columns=col)

    for gene in df['gene'].unique():
        geneDF = df[df['gene'] == gene]

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
                    intronDF = pd.concat([intronDF, row[col].to_frame().transpose()], ignore_index=True)
                    if count == -1: #start of table
                        count = 1
                        temp = pd.DataFrame(columns=col)
                        temp.loc[count,'gene'] = row['gene']
                        temp.loc[count,'chr'] = row['chr']
                        temp.loc[count,'stop'] = row['stop']

                        position = neg.loc[index, 'qstart']
                        prev_stop = row['stop']
                        check = 1
                        neg.loc[index,'check'] = check
                    elif index == len(neg)-1: #new gene end of table
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        check = 1
                        neg.loc[index,'check'] = check
                        neg = neg.drop([index])
                    else: # new gene
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        temp = pd.DataFrame(columns=col)
                        count = 1
                        temp.loc[count,'gene'] = row['gene']
                        temp.loc[count,'chr'] = row['chr']
                        temp.loc[count,'stop'] = row['stop']
                        position = neg.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        neg.loc[index,'check'] = check
                elif row['start']-prev_stop > 1000000: #large gap
                    intronDF = pd.concat([intronDF, row[col].to_frame().transpose()], ignore_index=True)
                    if index == len(neg)-1: #new gene end of table
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        check = 1
                        neg.loc[index,'check'] = check
                        neg = neg.drop([index])
                    else: #new gene
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        temp = pd.DataFrame(columns=col)
                        count += 1
                        temp.loc[count,'gene'] = row['gene']
                        temp.loc[count,'chr'] = row['chr']
                        temp.loc[count,'stop'] = row['stop']
                        position = neg.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        neg.loc[index,'check'] = check
                else: #next exon
                    intronDF = pd.concat([intronDF, row[col].to_frame().transpose()], ignore_index=True)
                    temp.loc[count,'start'] = row['start']
                    position = neg.loc[index, 'qstart']
                    prev_stop = row['stop']

                    check += 1
                    neg.loc[index,'check'] = check
                    if index == len(neg)-1: #end of table
                        temp.loc[count,'direction'] = '-'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

            intronDF = pd.concat([intronDF,neg])

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
                    intronDF = pd.concat([intronDF, row[col].to_frame().transpose()], ignore_index=True)
                    if count == -1: #start of table
                        count = 1
                        temp = pd.DataFrame(columns=col)
                        temp.loc[count,'gene'] = row['gene']
                        temp.loc[count,'chr'] = row['chr']
                        temp.loc[count,'start'] = row['start']

                        position = pos.loc[index, 'qstart']
                        prev_stop = row['stop']
                        check = 1
                        pos.loc[index,'check'] = check
                    elif index == len(pos)-1: #new gene end of table
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        check = 1
                        pos.loc[index,'check'] = check
                        pos = pos.drop([index])
                    else: # new gene
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        temp = pd.DataFrame(columns=col)
                        count = 1
                        temp.loc[count,'gene'] = row['gene']
                        temp.loc[count,'chr'] = row['chr']
                        temp.loc[count,'start'] = row['start']
                        position = pos.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        pos.loc[index,'check'] = check
                elif row['start']-prev_stop > 1000000: #large gap
                    intronDF = pd.concat([intronDF, row[col].to_frame().transpose()], ignore_index=True)
                    if index == len(pos)-1: #new gene end of table
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        check = 1
                        pos.loc[index,'check'] = check
                        pos = pos.drop([index])
                    else: #new gene
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])

                        temp = pd.DataFrame(columns=col)
                        count +=1
                        temp.loc[count,'gene'] = row['gene']
                        temp.loc[count,'chr'] = row['chr']
                        temp.loc[count,'start'] = row['start']
                        position = pos.loc[index, 'qstart']
                        prev_stop = row['stop']

                        check = 1
                        pos.loc[index,'check'] = check
                else: #next exon
                    intronDF = pd.concat([intronDF, row[col].to_frame().transpose()], ignore_index=True)
                    temp.loc[count,'stop'] = row['stop']
                    position = pos.loc[index, 'qstart']
                    prev_stop = row['stop']

                    check += 1
                    pos.loc[index,'check'] = check
                    if index == len(pos)-1: #end of table
                        temp.loc[count,'direction'] = '+'
                        temp.loc[count,'length'] = temp.loc[count,'stop'] - temp.loc[count,'start']
                        exonDF = pd.concat([exonDF,temp])


    intronDF = intronDF.sort_values(by=['start'])
    intronDF = intronDF.reset_index(drop=True)
    exonDF = exonDF.sort_values(by=['start'])
    exonDF = exonDF.reset_index(drop=True)
    return (exonDF,intronDF)

def drawGeneOutline(df,regionStart,regionEnd,scaleFactor,fig):
    # Reformat and store gene coordinates
    coordDict = {}
    regionStart = df['start'].min()
    regionEnd = df['stop'].max()
    regionLen = regionEnd - regionStart
    scaledRegionLen = regionLen / scaleFactor

    strandInfo = {}
    for index, row in df.iterrows():
        strandInfo[df.loc[index, 'gene']] = df.loc[index, 'direction']

        geneID = row['gene']
        geneDesc = row['stop']
        scaffoldID = row['chr']

        geneStart = row['start']
        geneEnd = row['stop']
        geneLen = geneEnd - geneStart

        newGeneStart = geneStart - regionStart + 1

        scaledGeneStart = newGeneStart / scaleFactor
        scaledGeneLen = geneLen / scaleFactor

        if (scaffoldID, regionLen, scaledRegionLen) not in coordDict:
            coordDict[(scaffoldID, regionLen, scaledRegionLen)] = []
        coordDict[(scaffoldID, regionLen, scaledRegionLen)].append((geneID, geneDesc, scaledGeneStart, scaledGeneLen))

    ### Draw outline  
    for scaffoldID, regionLen, scaledRegionLen in coordDict:
        for blastID, blastDesc, scaledGeneStart, scaledGeneLen in coordDict[(scaffoldID, regionLen, scaledRegionLen)]:
            x = scaledGeneStart+15
            y = 37

            if strandInfo[blastID] == '+':
                points = [(x-0.7, y),(x+scaledGeneLen+0.7, y),(x+scaledGeneLen+7, y+25.5),(x+scaledGeneLen+0.7, y+51),(x-0.7, y+51)]
            elif strandInfo[blastID] == '-':
                points = [(x-5.7, y+25.5),(x-0.7, y),(x+scaledGeneLen+0.7, y),(x+scaledGeneLen+0.7, y+51),(x-0.7, y+51)]
            fig.add(fig.polygon(points, fill='#8C8C8C', stroke='black'))

    return(fig)

def drawCassette(df, scaleFactor, name, chrom, shape, pad, dfI):
    # Reformat and store gene coordinates
    coordDict = {}
    regionStart = df['start'].min()
    regionEnd = df['stop'].max()
    regionLen = regionEnd - regionStart
    scaledRegionLen = regionLen / scaleFactor

    strandInfo = {}
    for index,row in df.iterrows():
        strandInfo[df.loc[index, 'gene']] = df.loc[index, 'direction']
        geneID = row['gene']
        geneDesc = row['stop']
        scaffoldID = row['chr']

        geneStart = row['start']
        geneEnd = row['stop']
        geneLen = geneEnd - geneStart

        newGeneStart = geneStart - regionStart + 1

        scaledGeneStart = newGeneStart / scaleFactor
        scaledGeneLen = geneLen / scaleFactor

        if (scaffoldID, regionLen, scaledRegionLen) not in coordDict:
            coordDict[(scaffoldID, regionLen, scaledRegionLen)] = []
        coordDict[(scaffoldID, regionLen, scaledRegionLen)].append((geneID, geneDesc, scaledGeneStart, scaledGeneLen))

    # Draw cassette
    outfile = name + "." + chrom + ".svg"
    figWidth = scaledRegionLen + 30
    figHeight = 250
    fig = svg.Drawing(
        filename=outfile,
        size=(figWidth,figHeight),
        profile='full'
    )

    for scaffoldID, regionLen, scaledRegionLen in coordDict:
        fig.add(fig.rect((0, 50), (figWidth, 25), fill='#d2d4d2', rx=2, ry=2))
        fig.add(fig.rect((0, 100), (30, 12), fill='#d2d4d2', rx=1, ry=1))
        desc = str(name) + "; " + chrom
        fig.add(fig.text(desc, insert=(40, 111), fill='black', font_size='16px'))

        if shape == 'pentaI':
            fig = drawGeneOutline(dfI,regionStart,regionEnd,scaleFactor,fig)

        for geneID, geneDesc, scaledGeneStart, scaledGeneLen in coordDict[(scaffoldID, regionLen, scaledRegionLen)]:
            x = scaledGeneStart + 15
            y = 37.5

            match shape:
                case 'rect':
                    points = [(x-pad,y), (x+scaledGeneLen+pad,y), (x+scaledGeneLen+pad,y+50), (x-pad,y+50)]

                case 'tri':
                    if strandInfo[geneID] == '+':
                        points = [(x, y), (x+scaledGeneLen+8, y+25), (x, y+50)]
                    elif strandInfo[geneID] == '-':
                        points = [(x, y), (x-scaledGeneLen-8, y+25), (x, y+50)]
                
                case 'trap':
                    if strandInfo[geneID] == '+':
                        points = [(x, y), (x+scaledGeneLen+8, y+10), (x+scaledGeneLen+8, y+40), (x, y+50)]
                    elif strandInfo[geneID] == '-':
                        points = [(x, y), (x-scaledGeneLen-8, y+10), (x-scaledGeneLen-8, y+40), (x, y+50)]

                case 'penta':
                    if strandInfo[geneID] == '+':
                        points = [(x-0.7, y),(x+scaledGeneLen+0.7, y),(x+scaledGeneLen+5.7, y+25.5),(x+scaledGeneLen+0.7, y+51),(x-0.7, y+51)]
                    elif strandInfo[geneID] == '-':
                        points = [(x-5.7, y+25.5),(x-0.7, y),(x+scaledGeneLen+0.7, y),(x+scaledGeneLen+0.7, y+51),(x-0.7, y+51)]

                case 'pentaI':
                    points = [(x,y), (x+scaledGeneLen,y), (x+scaledGeneLen,y+50), (x,y+50)]

                case 'flag':
                    y = 63
                    if strandInfo[geneID] == '+':
                        points = [(x,y),(x+scaledGeneLen+pad,y-25),(x-scaledGeneLen-pad,y-25)]
                    elif strandInfo[geneID] == '-':
                        points = [(x,y),(x+scaledGeneLen+pad,y+25),(x-scaledGeneLen-pad,y+25)]

            if shape == 'pentaI':
                fig.add(fig.polygon(points, fill='#e8000b'))
            else:
                fig.add(fig.polygon(points, fill='#e8000b', stroke='black'))

    fig.save()
    return()

def main():
    args = parse_arguments()

    match args.type:
        case 'bed':
            df = pd.read_csv(
                args.file,
                sep='\t',
                names=['chr','start','stop','gene','score','direction']
            )
            
            chrs = list(df['chr'].unique())
            for c in chrs:
                temp = df[df['chr']==c]
                temp = temp.sort_values(by=['start'])
                drawCassette(temp,args.scale,args.name,c,args.shape,args.pad,None)

        case 'out':
            df = pd.read_csv(
                args.file,
                sep='\t',
                header=None,
                names=['gene','chr', 'identity','length','mismatch','gaps','qstart','qstop','tstart','tstop','eval','bitscore']
            )
            
            chrs = list(df['chr'].unique())
            for c in chrs:
                temp = df[(df['chr'] == c)]
                if len(temp) > 10:
                    temp = getDIR(temp)
                    temp = temp.sort_values(by=['start'])
                    temp = temp.reset_index(drop=True)
                    
                    exonDF,intronDF = cleanDF(temp)
                    if not exonDF.empty:
                        if args.shape == 'pentaI':
                            drawCassette(intronDF,args.scale,args.name,c,args.shape,args.pad,exonDF)
                        else:
                            drawCassette(exonDF,args.scale,args.name,c,args.shape,args.pad,intronDF)

        case 'gff':
            df = pd.read_csv(
                args.file,
                sep='\t',
                header=None,
                names=['chr','source','type','start','stop','score','direction','phase','attribute']
            )
            df['gene'] = df['attribute'].str.split(";").str[0].str.split("=").str[1].str.replace('.cds','')
            df = df.drop(columns=['source','score','phase','attribute'])

            chrs = list(df['chr'].unique())
            for c in chrs:
                temp = df[df['chr']==c]
                temp = temp.sort_values(by=['start'])

                exonDF = temp[temp['type']=='CDS']
                intronDF = temp[temp['type']=='mRNA']
                if not exonDF.empty:
                    if args.shape == 'pentaI':
                        drawCassette(exonDF,args.scale,args.name,c,args.shape,args.pad,intronDF)
                    else:
                        drawCassette(intronDF,args.scale,args.name,c,args.shape,args.pad,exonDF)




if __name__ == '__main__':
    main()