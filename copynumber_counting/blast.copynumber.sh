#!/bin/bash

q_array=( "AAE1.aa.fasta" "OAC.aa.fasta" "PT4.aa.fasta" "OLS.aa.fasta" "GPPS_ls.aa.fasta" "GPPS_ss.aa.fasta" "BKR.cds.fasta" "ALT4.cds.fasta" )

FILE=$1
NAME=${FILE%%.*}
TYPE=$2

if [ ! -f ${NAME}.softmasked.fasta ]; then
	aws s3 cp s3://salk-tm-shared/csat/releases/${TYPE}/${NAME}/${FILE} .
	aws s3 cp s3://salk-tm-dev/allen_temp/cannabis/synthase_cassettes/csv_output/${NAME}_filterhits.csv .
	if [ "${FILE##*.}" = "gz" ]; then
		echo "##### UNZIPPING FILE"
		gunzip $FILE
		FILE=${NAME}.softmasked.fasta
	fi
else
	FILE=${NAME}.softmasked.fasta
fi

makeblastdb -in ${FILE} -dbtype nucl

for QFILE in "${q_array[@]}"
do
	QUERY=${QFILE%%.*}
	if [[ "${QFILE}" == *'aa'* ]]; then
		tblastn -query ${QFILE} -db ${FILE} -out ${NAME}.${QUERY}.tblastn.out -outfmt 6
	elif [[ "${QFILE}" == *'cds'* ]]; then
		blastn -query ${QFILE} -db ${FILE} -out ${NAME}.${QUERY}.blastn.out -outfmt 6
	fi
done

python count.copynumber.py ${NAME}

rm ${NAME}.*
rm ${NAME}_filterhits.csv
