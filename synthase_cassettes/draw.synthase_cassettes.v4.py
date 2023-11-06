import pandas as pd
import svgwrite as svg
import sys

##### FILTER SYNTHASE HITS #####
def filterSYN(syn):
    # Filter params
    bitscore = 0
    length = 500

    # Get synthase hits
    sdf = pd.read_csv(syn, sep='\t', header=None)
    sdf.columns = ['query', 'target', 'pident', 'length', 'mismatch', 'gaps', \
                  'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore']
    sdf = sdf[(sdf['bitscore'] >= bitscore) & (sdf['length'] >= length)] #filter by params

    # Add type and direction
    sdf['type'] = sdf[['mismatch', 'gaps']].apply(lambda row: 'Full' if (row['mismatch']<10 and row['gaps']==0) else 'Partial', axis=1) 
    for index, row in sdf.iterrows():
        if row['tstart'] < row['tend']:
            sdf.loc[index, 'direction'] = '+'
        else:
            sdf.loc[index, 'direction'] = '-'
    
    return(sdf)

##### FILTER LTR HITS #####
def filterLTR(ltr, targets):
    # Parameters
    bitscore = 1250

    # Start empty table
    ldf = pd.DataFrame(columns = ['query', 'target', 'pident', 'length', 'mismatch', 'gaps', \
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore'])

    # Concat rows with target chroms
    df = pd.read_csv(ltr, sep='\t', header=None)
    df.columns = ['query', 'target', 'pident', 'length', 'mismatch', 'gaps', \
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore']
    
    for target in targets:
        tempdf = df[(df['target'] == target)  & (df['bitscore'] >= bitscore) & (df['query'] == 'LTR08_09')] #filter by param
        ldf = pd.concat([ldf, tempdf])
    
    # Add type and direction
    ldf['type'] = 'LTR'
    for index, row in ldf.iterrows():
        if row['tstart'] < row['tend']:
                ldf.loc[index, 'direction'] = '+'
        else:
                ldf.loc[index, 'direction'] = '-'

    return(ldf)

##### RENAME HITS #####
def renameHITS(sdf, ldf, targets):
    # Max distance of CBDAS from LTR08
    dist = 60000
    
    fdf = pd.DataFrame(columns = ['query', 'target', 'pident', 'length', 'mismatch', 'gaps', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore', 'start', 'end', 'down', 'up'])
    fd = {}
    
    # Loop over target chromosomes
    for target in targets:
        tempSDF = sdf[sdf['target'] == target]
        
        # Start an empty table
        df = pd.DataFrame(columns = ['query', 'target', 'pident', 'length', 'mismatch', 'gaps', \
            'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore', 'start', 'end', 'down', 'up'])
        
        # Loop over synthase hits and rename correctly
        for index, row in tempSDF.iterrows():
            tempLDF = ldf.copy() #copy ltr table
            tempLDF.loc[len(tempLDF.index)] = row #add synthase to table

            # Sort and reorder coords
            tempLDF.sort_values('tstart', inplace=True)
            for index, row in tempLDF.iterrows():
                if row['direction'] == '-':
                    tempLDF.loc[index, 'start'] = row['tend']
                    tempLDF.loc[index, 'end'] = row['tstart']
                else:
                    tempLDF.loc[index, 'start'] = row['tstart']
                    tempLDF.loc[index, 'end'] = row['tend']

            # Calculate distance
            tempLDF['down'] = tempLDF['start'] - tempLDF['end'].shift(1)
            tempLDF['up'] = tempLDF['start'].shift(-1) - tempLDF['end']
            tempLDF = tempLDF[tempLDF['query'] != 'LTR08_09']
            fdf = pd.concat([fdf, tempLDF])

            # Rename hit
            for index, row in tempLDF.iterrows():
                if row['type'] == 'Full' and row['query'] == 'Cannabis_sativa_CBDAS_AB292682':
                    tempLDF.loc[index,'query'] = "CBDAS"
                elif row['type'] == 'Full' and row['query'] == 'Cannabis_sativa_THCAS_AB057805':
                    tempLDF.loc[index,'query'] = "THCAS"
                elif row['type'] == 'Full' and row['query'] == "Cannabis_sativa_CBCAS_LY658671.1":
                    tempLDF.loc[index,'query'] = "CBCAS"
                elif row['down'] < dist or row['up'] < dist:
                    tempLDF.loc[index,'query'] = "CBDAS"
            df = pd.concat([df ,tempLDF])
             
        # Sort new table
        df.sort_values('tstart',inplace=True)
        df.reset_index(drop=True, inplace=True)
        fd[target] = removeDUPES(df)
    return(fd, fdf)

##### REMOVE DUPLICATES #####
def removeDUPES(df):
    # Get hits with the same starting position
    startDict = {}
    for index, row in df.iterrows():
        if row['start'] in startDict:    
            startDict[row['start']].append(index)
        else:
            startDict[row['start']] = []
            startDict[row['start']].append(index)
    # Remove duplicate hits
    for key in startDict.keys():
        if len(startDict[key]) < 2:
            pass
        else:
            keep = "none"
            for index in startDict[key]:
                if keep == "none":
                    keep = index
                elif df.loc[index, 'type'] == "Full":
                    df.drop(keep, inplace=True)
                    keep = index                            
                elif df.loc[index, 'bitscore'] > df.loc[keep, 'bitscore']:
                    df.drop(keep, inplace=True)
                    keep = index
                else:
                    df.drop(index, inplace=True)

    # Get hits with the same ending position
    endDict = {}
    for index, row in df.iterrows():
        if row['end'] in endDict:    
            endDict[row['end']].append(index)
        else:
            endDict[row['end']] = []
            endDict[row['end']].append(index)
    # Remove duplicate hits
    for key in endDict.keys():
        if len(endDict[key]) < 2:
            pass
        else:
            keep = "none"
            for index in endDict[key]:
                if keep == "none":
                    keep = index
                elif df.loc[index, 'type'] == "Full":
                    df.drop(keep, inplace=True)
                    keep = index
                elif df.loc[index, 'bitscore'] > df.loc[keep, 'bitscore']:
                    df.drop(keep, inplace=True)
                    keep = index
                else:
                    df.drop(index, inplace=True)
                    
    # Rename THCAS and CBDAS
    for index, row in df.iterrows():
        if row['query'] == 'Cannabis_sativa_THCAS_AB057805':
            df.loc[index,'query'] = "THCAS"
        elif row['query'] == "Cannabis_sativa_CBCAS_LY658671.1":
            df.loc[index,'query'] = "CBCAS"
            
    return(df)

##### MAKE FIGURES #####
def makeFIG(fdf, target, count, name):
    # Coordinate to figure scaling
    scaleFactor = 2000
    
    # Rename hits
    df = fdf.reset_index(drop=True)
    for index, row in df.iterrows():
        df.loc[index, 'query'] = df.loc[index, 'query'] + "_" + str(index)
        
    # Save hit direction
    strandInfo = {}
    for index, row in df.iterrows():
        strandInfo[df.loc[index, 'query']] = df.loc[index, 'direction']
        
    # Reformat and save hit coordinates
    blastCoordDict = {}

    regionStart = df['start'].min() #region coords
    regionEnd = df['end'].max()
    regionLen = regionEnd - regionStart
    scaledRegionLen = regionLen / scaleFactor

    newRegionStart = 1 #adjusted region coords
    newRegionEnd = newRegionStart + regionLen

    for index, row in df.iterrows():
        blastID = row['query'] #blast hit
        blastDesc = row['type']
        scaffoldID = row['target']

        geneStart = row['start'] #hit coords
        geneEnd = row['end']
        geneLen = geneEnd - geneStart

        newGeneStart = geneStart - regionStart + 1 #adjusted hit coords
        newGeneEnd = newGeneStart + geneLen + 1

        scaledGeneStart = newGeneStart / scaleFactor #scale hits for smaller plot
        # scaledGeneEnd = newGeneEnd / scaleFactor
        scaledGeneLen = geneLen / scaleFactor

        if (scaffoldID, regionLen, scaledRegionLen) not in blastCoordDict:
            blastCoordDict[(scaffoldID, regionLen, scaledRegionLen)] = []
        blastCoordDict[(scaffoldID, regionLen, scaledRegionLen)].append((blastID, blastDesc, scaledGeneStart, scaledGeneLen))
    
    ### Draw figure
    # geneCount = df.shape[0]
    figWidth = scaledRegionLen+30
    figHeight = 250
    file = name + "_" + str(count) + ".svg"
    fig = svg.Drawing(filename=file, size=(figWidth, figHeight), profile='full')
        
    for scaffoldID, regionLen, scaledRegionLen in blastCoordDict:
        fig.add(fig.rect((0, 50), (figWidth, 25), fill='#d2d4d2', rx=2, ry=2))
        fig.add(fig.rect((0, 100), (30, 12), fill='#d2d4d2', rx=1, ry=1))
        desc = str(name) + "; " + target# + "; " + str(int(regionLen)) + "bp"
        fig.add(fig.text(desc, insert=(40, 111), fill='black', font_size='16px'))
        for blastID, blastDesc, scaledGeneStart, scaledGeneLen in blastCoordDict[(scaffoldID, regionLen, scaledRegionLen)]:
            x = scaledGeneStart + 10
            y = 37.5
            if 'THCAS' in blastID:
                if 'Full' not in blastDesc:
                    if strandInfo[blastID] == '+':
                        points = [(x, y), (x+scaledGeneLen+8, y+10), (x+scaledGeneLen+8, y+40), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#e164f0', stroke='black'))
                    else:
                        points = [(x, y), (x-scaledGeneLen-8, y+10), (x-scaledGeneLen-8, y+40), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#e164f0', stroke='black'))

                else:
                    if strandInfo[blastID] == '+':
                        points = [(x, y), (x+scaledGeneLen+15, y+25), (x, y+50)]    
                        fig.add(fig.polygon(points, fill='#b41964', stroke='black'))
                    else:
                        points = [(x, y), (x-scaledGeneLen-15, y+25), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#b41964', stroke='black'))
                            
            elif 'CBCAS' in blastID:
                if 'Full' not in blastDesc:
                    if strandInfo[blastID] == '+':
                        points = [(x, y), (x+scaledGeneLen+8, y+10), (x+scaledGeneLen+8, y+40), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#ffd700', stroke='black'))
                    else:
                        points = [(x, y), (x-scaledGeneLen-8, y+10), (x-scaledGeneLen-8, y+40), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#ffd700', stroke='black'))

                else:
                    if strandInfo[blastID] == '+':
                        points = [(x, y), (x+scaledGeneLen+15, y+25), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#e17800', stroke='black'))
                    else:
                        points = [(x, y), (x-scaledGeneLen-15, y+25), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#e17800', stroke='black'))

            else:
                if 'Full' not in blastDesc:
                    if strandInfo[blastID] == '+':
                        points = [(x, y), (x+scaledGeneLen+8, y+10), (x+scaledGeneLen+8, y+40), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#00c8f0', stroke='black'))
                    else:
                        points = [(x, y), (x-scaledGeneLen-8, y+10), (x-scaledGeneLen-8, y+40), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#00c8f0', stroke='black'))

                else:
                    if strandInfo[blastID] == '+':
                        points = [(x, y), (x+scaledGeneLen+15, y+25), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#0064e1', stroke='black'))
                    else:
                        points = [(x, y), (x-scaledGeneLen-15, y+25), (x, y+50)]
                        fig.add(fig.polygon(points, fill='#0064e1', stroke='black'))
    fig.save()
    return()


def main():
    syn = sys.argv[1]
    ltr = sys.argv[2]
    name = syn.split('.')[0]

    sdf = filterSYN(syn)
    targets = list(sdf['target'].unique())
    
    ldf = filterLTR(ltr, targets)

    fd, fdf = renameHITS(sdf, ldf, targets)

    pdf = pd.DataFrame(columns = ['query', 'target', 'pident', 'length', 'mismatch', 'gaps', \
                        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bitscore', 'start', 'end', 'down', 'up'])

    for target in targets:
        if fd[target].shape[0] < 4:
            print("SKIPPING " + str(target) + " : " + str(fd[target].shape[0]))
        else:
            makeFIG(fd[target], target, list(fd.keys()).index(target), name)
            pdf = pd.concat([pdf, fd[target]], ignore_index=True)

    fdf.to_csv(str(name) + '_allhits.csv')
    pdf.to_csv(str(name) + '_filterhits.csv')

if __name__ == '__main__':
    main()