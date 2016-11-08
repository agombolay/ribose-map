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
	samtools idxstats $input | cut -f 1,3 > $output
	
	#Select only chromosomes of interest (Remove reads that align to unidentified regions of genome)
	#grep -v 'chr[a-zA-Z0-9]\+_[a-zA-Z0-9]\+' temporary.txt | grep -v '*' | grep -v 'chrEBV' > $output
	
	#Calculate total number of reads in genome
	#cat $output | awk '{sum+=$2} END{print "Total =",sum}'
	
	#Calculate number of reads that align to chromosomes
	#grep -v 'chrM' $output | awk '{sum+=$2} END{print "Nuclear =",sum}'
	
	#Calculate number of reads that align to mitochondria
	#grep 'chrM' $output | awk '{sum+=$2} END{print "Mitochondria =",sum}'
	
done

#Remove temporary file
rm temporary.txt
