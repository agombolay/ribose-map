#!/usr/bin/env bash

#Author: Alli Gombolay
#This program converts the raw sequencing files from SRA (.sra) format to FASTQ (.fq)

#COMMAND LINE OPTIONS

#Name of the program (1_setUpRawData.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-i] '/path/to/file1.sra etc.' [-h]
          -i Filepaths of input .sra sequencing files"
}

#Use getopts function to create the command-line options ([-i] and [-h])
while getopts "i:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) sra=($OPTARG) ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

#Convert raw sequencing files from .sra to .fq
for samples in ${sra[@]};
do
        #Extract input directories from filepaths
	inputDirectory=$(dirname "${sra}")
	
	#Change current diretory to 
	cd $inputDirectory
	
        fastq-dump $samples
done
