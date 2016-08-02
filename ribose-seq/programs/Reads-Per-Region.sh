#! /usr/bin/env bash

#Author: Alli Gombolay
#Date: August 2, 2016
#This program extracts number of reads per chromosome and saves output to text file in same location as bam file

#COMMAND LINE OPTIONS

#Name of the program (Reads-Per-Region.sh)
program=$0

#Usage statement of the program
function usage () {
	echo "Usage: $program [-i] '/path/to/file1.bam etc.' [-d] 'Ribose-seq directory' [-h]
	-i Input (i.e., /projects/home/agombolay3/data/repository/Ribose-seq-Project/ribose-seq/results/hg38/FS56/Alignment)
	-d Location of user's local Ribose-seq directory (i.e., /projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Use getopts function to create the command-line options ([-i], [-d], and [-h])
while getopts "i:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) bam=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

for samples in ${bam[@]};
do
	#Extract sample names from filepaths
	filename=$(basename "${samples}")
	samples="${filename%.*}"
	
	#Extract input directories from filepaths
	inputDirectory=$(dirname "${bam}")
	
	#INPUT
	#Location of input files
	input=$inputDirectory/$samples.bam

	#OUTPUT
	#Location of output files
	output=$inputDirectory/$samples.Reads-Per-Region.txt

	samtools idxstats $samples.bam | cut -f 1,3 > $output
done
