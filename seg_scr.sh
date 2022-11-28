#!/bin/bash

INPUT_DIRECTORY=$1
OUTPUT_DIRECTORY=$2
PROGRAM=seg

if [ ! -d $OUTPUT_DIRECTORY ]; then
	mkdir $OUTPUT_DIRECTORY
fi

for FILE in $INPUT_DIRECTORY/*.fasta; do 
	FILENAME=$(echo "$(basename $FILE)" | cut -f 1 -d '.')
	echo Making a file: ${OUTPUT_DIRECTORY}/$FILENAME
	$PROGRAM $FILE > "${OUTPUT_DIRECTORY}/${FILENAME}.txt"; 
done

