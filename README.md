# cannabis_pangenome
```
conda create -c conda-forge -c bioconda -c anaconda -n cassettes "python>=3.10" pandas svgwrite blast seaborn bedtools mafft fasttree awscli
```

### synthase_cassettes
Directory containing files to generate gene cassettes for the cannabinoid synthases (THCAS, CBDAS, CBCAS).
Output svg and csv files will be uploaded to s3://salk-tm-dev/allen_temp/synthase_cassettes/
```
cd synthase_cassettes

nohup sh run.draw_synthase_cassettes.sh publics.genomes.list &> publics.out &

nohup sh run.draw_synthase_cassettes.sh scaffolded.genomes.list &> scaffolded.out &

nohup sh run.draw_synthase_cassettes.sh not_scaffolded.genomes.list &> not_scaffolded.out &
```

### copynumber_counting
Directory containing files to get copy numbers for the pathway genes (AAE1, OAC, PT4, OLS) and synthases.
Counts will be appended to 'pathway_copynumbers.csv'.
```
cd copynumber_counting

nohup sh run.count_copynumbers.sh publics.genomes.list &> publics.out &

nohup sh run.count_copynumbers.sh scaffolded.genomes.list &> scaffolded.out &

nohup sh run.count_copynumbers.sh not_scaffolded.genomes.list &> not_scaffolded.out &

python draw.copynumber.py pathway_copynumbers.csv
```

### gene_tree
Directory containing files to obtain fasta files from blast tables and generate a gene tree.
Fasta files will be placed in the "fastas" directory. Gene tree will output as "full_synthases.tree".
```
cd gene_tree

mkdir fastas

nohup sh run.getfasta.sh publics.genomes.list &> publics.out &

nohup sh run.getfasta.sh scaffolded.genomes.list &> scaffolded.out &

nohup sh run.getfasta.sh not_scaffolded.genomes.list &> not_scaffolded.out &

cat fastas/*.fasta > full_synthases.fasta

mafft --auto --reorder full_synthases.fasta > full_synthases.aln

FastTree -nt full_synthases.aln > full_synthases.tree
```