#!/bin/bash

for path in AH3Mb_alignments/*.paf; do
    filename=$(basename "$path")
    genome="${filename#*.*.}"
    genome="${genome%%.*}"

    grep "AH3Mb.chrY" "$path" > temp.txt
    grep "${genome}.chrY" temp.txt > "AH3Mb.${genome}.chrY.paf"
done

for path in AH3Ma_alignments/*.paf; do
    filename=$(basename "$path")
    genome="${filename#*.*.}"
    genome="${genome%%.*}"

    grep "AH3Ma.chrX" "$path" > temp.txt
    grep "${genome}.chrX" temp.txt > "AH3Ma.${genome}.chrX.paf"
done

rm temp.txt

cat AH3Mb.*.chrY.paf > chrY.paf
cat AH3Ma.*.chrX.paf > chrX.paf

rm AH3Mb.*
rm AH3Ma.*