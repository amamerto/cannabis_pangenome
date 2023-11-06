#!/bin/bash

LIST=$1
TEMP=${LIST##*/}
TYPE=${TEMP%%.*}

while read FILE; do
    sh getfasta.from_csv.sh ${FILE} ${TYPE}
done < ${LIST}
