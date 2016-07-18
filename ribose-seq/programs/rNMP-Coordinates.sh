#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: July 11, 2016
#This program determines 0-based and 1-based coordinates of rNMPs (5’ end of each mapped read) for positive and negative strands

#COMMAND LINE OPTIONS

#Name of the program (rNMP-Coordinates.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-i] '/path/to/file1.bam etc.' [-d] 'Ribose-seq directory' [-h]
          -i Filepaths of input BAM files
          -d Location of user's local Ribose-seq directory"
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

#Obtain 0-based and 1-based coordinates of rNMPs
for samples in ${bam[@]};
do

	#Extract sample names from filepaths
	filename=$(basename "${bam}")
	samples="${filename%.*}"

	#Extract input directories from filepaths
	inputDirectory=$(dirname "${bam}")

	#INPUT
	#Location of input BAM files
	input=$inputDirectory/$samples.bam

	#OUTPUT
	#Location of output "ribose-seq" alignment directory
	output=$directory/ribose-seq/results/align_sacCer2/$samples/alignment

	bed=$output/$samples.bed
	sam=$output/$samples.sam
	fasta=$output/$samples.fasta

	#COORDINATES (0-BASED) of SEQUENCING READS

	#Covert file from BAM to BED format using BEDtools software
	bedtools bamtobed -i $input > $bed

	#Convert file from BAM to SAM format using SAMtools software
	samtools view $input > $sam

	#Extract read coordinates, sequences, and strands from BED and SAM files and save it to new file
	paste $bed $sam $fasta | awk -v "OFS=\t" '{print $1, $2, $3, $4, $16, $6, $21}' > $output/$samples.coordinates.bed

	#0-BASED COORDINATES OF rNMPs

	#Obtain positions of rNMPs (3’ end of each mapped read) for positive strand:
	bedtools genomecov -3 -strand + -bg -ibam $input > $output/$samples.rNMPs.positive.0-based.txt

	#Obtain positions of rNMPs (3’ end of each mapped read) for negative strand:
	bedtools genomecov -3 -strand - -bg -ibam $input > $output/$samples.rNMPs.negative.0-based.txt

	#1-BASED COORDINATES OF	rNMPs

	#Obtain positions of rNMPs (3’ end of each mapped read) for positive strand:
	bedtools genomecov -3 -strand + -d -ibam $input > $output/$samples.rNMPs.positive.1-based.txt

	#Remove rows where genome coverage equals 0
	awk '$3 != 0' $output/$samples.rNMPs.positive.1-based.txt > temporary

	#Change filename back to original
	mv temporary $output/$samples.rNMPs.positive.1-based.txt

	#Obtain positions of rNMPs (3’ end of each mapped read) for negative strand:
	bedtools genomecov -3 -strand - -d -ibam $input > $output/$samples.rNMPs.negative.1-based.txt

	#Remove rows where genome coverage equals 0
	awk '$3 != 0' $output/$samples.rNMPs.negative.1-based.txt > temporary

	#Change filename back to original
	mv temporary $output/$samples.rNMPs.negative.1-based.txt

done


