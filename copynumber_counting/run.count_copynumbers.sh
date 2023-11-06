#!/bin/bash

LIST=$1
TEMP=${LIST##*/}
TYPE=${TEMP%%.*}

while read FILE; do
    bash blast.copynumber.sh ${FILE} ${TYPE}
done < ${LIST}
