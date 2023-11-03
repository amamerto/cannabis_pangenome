#!/bin/bash

FILE=$1.softmasked.fasta.gz
NAME=${FILE%%.*}
TYPE=$2

echo "~~~#####   STARTING ${NAME}   #####~~~"
aws s3 cp s3://salk-tm-shared/csat/releases/$TYPE/$NAME/$FILE .
if [ "${FILE##*.}" = "gz" ]; then
        echo "##### UNZIPPING FILE"
        gunzip $FILE
        FILE=${NAME}.softmasked.fasta
fi

echo "~~~#  GRABBING ${NAME} FILTERED HITS  #~~~"
aws s3 cp s3://salk-tm-dev/allen_temp/synthase_cassettes/csv_output/${NAME}_filterhits.csv .

echo "~~~#####   MAKING ${NAME} BED   #####~~~"
python filter_blastn.py ${NAME}_filterhits.csv

echo "~~~#####   MAKING ${NAME} FASTA   #####~~~"
bedtools getfasta -name -fi ${NAME}.softmasked.fasta -bed ${NAME}_fullsyn.bed -fo ${NAME}_fullsyn.fasta
sed -i 's/\(.*\)*::\(.*\)/\1/g' ${NAME}_fullsyn.fasta
mv ${NAME}_fullsyn.fasta fastas/${NAME}_fullsyn.fasta

rm ${NAME}_fullsyn.bed
rm ${NAME}_filterhits.csv
rm ${NAME}.softmasked.fasta*