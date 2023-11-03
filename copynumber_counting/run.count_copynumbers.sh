#!/bin/bash

LIST=$1
TEMP=${LIST##*/}
TYPE=${TEMP%%.*}

while read FILE; do
    sh blast.copynumber.sh ${FILE} ${TYPE}
done < ${LIST}
