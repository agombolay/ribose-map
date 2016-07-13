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
for files in ${bam[@]};
do

	#COORDINATES (0-BASED) of SEQUENCING READS

	#Covert file from BAM to BED format using BEDtools software
	bedtools bamtobed -i $files.bam > $files.bed

	#Convert file from BAM to SAM format using SAMtools software
	samtools view $files.bam > $files.sam

	#Extract read coordinates, sequences, and strands from BED and SAM files and save it to new file
	paste $files.bed $files.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' > $files.coordinates.bed

	#0-BASED COORDINATES OF rNMPs

	#Obtain positions of rNMPs (5’ end of each mapped read) for positive strand:
	bedtools genomecov -5 -strand + -bg -ibam $files.bam > $files.rNMPs.positive.txt

	#Obtain positions of rNMPs (5’ end of each mapped read) for negative strand:
	bedtools genomecov -5 -strand - -bg -ibam $files.bam > $files.rNMPs.negative.txt

	#1-BASED COORDINATES OF	rNMPs

	#Obtain positions of rNMPs (5’ end of each mapped read) for positive strand:
	bedtools genomecov -5 -strand + -d -ibam $files.bam > $files.rNMPs.positive.txt

	#Remove rows where genome coverage equals 0
	awk '$3 != 0' $files.rNMPs.positive.txt > temporary

	#Change filename back to original
	mv temporary $files.rNMPs.positive.txt

	#Obtain positions of rNMPs (5’ end of each mapped read) for negative strand:
	bedtools genomecov -5 -strand - -d -ibam $files.bam > $files.rNMPs.negative.txt

	#Remove rows where genome coverage equals 0
	awk '$3 != 0' $files.rNMPs.negative.txt > temporary

	#Change filename back to original
	mv temporary $files.rNMPs.negative.txt

done


