#!/bin/bash

chrom=$1

for path in other_pafs/*.paf; do
    filename=$(basename "$path")
    genome="${filename#*.*.}"
    genome="${genome%%.*}"

    grep "AH3Mb.chrY" "$path" > temp.txt
    grep "${genome}.${chrom}" temp.txt > "AH3Mb.${genome}.${chrom}.paf"
done

rm temp.txt

cat AH3Mb.*.${chrom}.paf > ${chrom}.paf