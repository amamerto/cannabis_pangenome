#!/bin/bash

LIST=$1
TEMP=${LIST##*/}
TYPE=${TEMP%%.*}

while read FILE; do
    sh get.fasta_fromcsv.sh ${FILE} ${TYPE}
done < ${LIST}
