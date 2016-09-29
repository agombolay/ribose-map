#! /usr/bin/env bash

#Author: Alli Gombolay
#Date: August 2, 2016
#This program extracts number of reads per chromosome and saves output to text file in same location as input bam file

#COMMAND LINE OPTIONS

#Name of the program (9_Reads-Per-Region.sh)
#program=$0

#Usage statement of the program
function usage () {
	echo "Usage: Reads-Per-Region.sh [-i] 'BAM' [-h]
		-i Filepath of BAM files ('/path/to/FS1.bam' etc.)"
}

#Use getopts function to create the command-line options ([-i] and [-h])
while getopts "i:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) bam=($OPTARG) ;;
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

	#Extract number of reads per chromosome and save output
	samtools idxstats $input | cut -f 1,3

done
