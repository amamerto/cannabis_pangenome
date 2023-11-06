#!/bin/bash

LIST=$1
TEMP=${LIST##*/}
TYPE=${TEMP%%.*}

while read FILE; do
    sh blast.synthase_cassettes.sh ${FILE} ${TYPE}
done < ${LIST}
