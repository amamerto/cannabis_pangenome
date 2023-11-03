#!/bin/bash

LIST=$1
TEMP=${LIST##*/}
TYPE=${TEMP%%.*}

while read FILE; do
    sh draw.synthase_cassettes.v4.sh ${FILE} ${TYPE}
done < ${LIST}
