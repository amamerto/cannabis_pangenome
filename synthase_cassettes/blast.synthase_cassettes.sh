#!/bin/bash

FILE=$1
NAME=${FILE%%.*}
TYPE=$2

# echo "~~~#####   DOWNLOADING ${NAME} ASM   #####~~~"
# aws s3 cp s3://salk-tm-shared/csat/releases/$TYPE/$NAME/$FILE .
# if [ "${FILE##*.}" = "gz" ]; then
#         echo "##### UNZIPPING FILE"
#         gunzip $FILE
#         FILE=${NAME}.softmasked.fasta
# fi

# echo "##### MAKING BLAST DB"
# makeblastdb -in $FILE -dbtype nucl
# echo "##### FINDING SYNTHASE HITS"
# blastn -query query.cds_synthases.fasta -db $FILE -out ${NAME}.blastn.syn -outfmt 6
# echo "##### FINDING LTR08 HITS"
# blastn -query query.ltr08.fasta -db $FILE -out ${NAME}.blastn.ltr -outfmt 6

echo "##### DRAWING CASSETTES"
python draw.synthase_cassettes.v4.py ${NAME}.blastn.syn ${NAME}.blastn.ltr

# echo "##### UPLOADING TO AWS"
# for SVG in *.svg; do
#         aws s3 cp $SVG s3://salk-tm-dev/allen_temp/synthase_cassettes/svg_output/
# done
# for CSV in *.csv; do
#         aws s3 cp $CSV s3://salk-tm-dev/allen_temp/synthase_cassettes/csv_output/
# done

# echo "##### REMOVING FILES"
# rm ${NAME}.blastn.*
# rm ${NAME}_*.svg
# rm ${NAME}.softmasked.fasta*
# rm ${NAME}_*.csv

# echo "~~~#####   COMPLETED ${NAME}   #####~~~"