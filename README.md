# cannabis_pangenome
```
conda create -c conda-forge -c bioconda -c anaconda -n cassettes "python>=3.10" pandas plotnine svgwrite blast seaborn bedtools mafft fasttree awscli
```

### synthase_cassettes
Directory containing files to generate gene cassettes for the cannabinoid synthases (THCAS, CBDAS, CBCAS).
Output svg and csv files will be uploaded to s3://salk-tm-dev/allen_temp/synthase_cassettes/
```
cd synthase_cassettes

nohup sh run.draw_synthase_cassettes.sh ../publics.genomes.list &> publics.out &

nohup sh run.draw_synthase_cassettes.sh ../scaffolded.genomes.list &> scaffolded.out &

nohup sh run.draw_synthase_cassettes.sh ../not_scaffolded.genomes.list &> not_scaffolded.out &
```

### copynumber_counting
Directory containing files to get copy numbers for the pathway genes (AAE1, OAC, PT4, OLS) and synthases.
Counts will be appended to 'pathway_copynumbers.csv'.
```
cd copynumber_counting

nohup sh run.count_copynumbers.sh ../publics.genomes.list &> publics.out &

nohup sh run.count_copynumbers.sh ../scaffolded.genomes.list &> scaffolded.out &

nohup sh run.count_copynumbers.sh ../not_scaffolded.genomes.list &> not_scaffolded.out &

python draw.copynumber.py pathway_copynumbers.csv
```

### gene_tree
Directory containing files to obtain fasta files from blast tables and generate a gene tree.
Fasta files will be placed in the "fastas" directory. Gene tree will output as "full_synthases.tree".
```
cd gene_tree

mkdir fastas

nohup sh run.getfasta.sh ../publics.genomes.list &> publics.out &

nohup sh run.getfasta.sh ../scaffolded.genomes.list &> scaffolded.out &

nohup sh run.getfasta.sh ../not_scaffolded.genomes.list &> not_scaffolded.out &

cat fastas/*.fasta > full_synthases.fasta

mafft --auto --reorder full_synthases.fasta > full_synthases.aln

FastTree -nt full_synthases.aln > full_synthases.tree
```

### pandots
Directory containing files for filtering single chromosomes from paf files and drawing multi-genome dotplots.
Ouput will be saved as svg.

Plot Autosome
```
cd pandots

#Run single chromosome
sh filter_autosomes.sh chr7

python pandots.py --paf chr7.paf --ref EH23a.chr7 --out chr7 --key True --reorder True --recolor True


#Or run all chromosome commands
for i in {1..9}; do
    sh filter_autosomes.sh chr${i}
    python pandots.py --paf chr${i}.paf --ref EH23a.chr${i} --out chr${i} --key True --reorder True --recolor True
done
```

Plot Allosome
```
cd pandots

sh filter_allosomes.sh

python pandots.py --paf chrX.paf --ref AH3Ma.chrX --out chrX --key True --invert True
python pandots.py --paf chrY.paf --ref AH3Mb.chrY --out chrY --key True --flip True
```