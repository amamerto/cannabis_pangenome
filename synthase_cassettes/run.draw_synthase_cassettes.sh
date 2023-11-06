#!/bin/bash

LIST=$1
TEMP=${LIST##*/}
TYPE=${TEMP%%.*}

while read FILE; do
    bash blast.synthase_cassettes.sh ${FILE} ${TYPE}
done < ${LIST}
