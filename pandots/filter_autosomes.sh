#!/bin/bash

chrom=$1

for path in EH23a_alignments/*.paf; do
    filename=$(basename "$path")
    genome="${filename#*.*.}"
    genome="${genome%%.*}"
    
    grep "EH23a.${chrom}" "$path" > temp.txt
    grep "${genome}.${chrom}" temp.txt > "EH23a.${genome}.${chrom}.paf"
done

rm temp.txt

cat EH23a.*.${chrom}.paf > ${chrom}.paf

rm EH23a.*