#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: July 11, 2016
#This program determines 0-based and 1-based coordinates of rNMPs (5’ end of each mapped read) for + and - strands

#COMMAND LINE OPTIONS

#Name of the program (rNMP-Coordinates.sh)
program=$0

#Usage statement of the program
function usage () {
	echo "Usage: Ribonucleotide-Coordinates.sh [-i] 'BAM' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Filepaths of BAM files ('/path/to/FS1.bam etc.')
	-r Name of reference genome folder in which to store output files (i.e., sacCer2)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
}

#Use getopts function to create the command-line options ([-i], [-r], [-d], and [-h])
while getopts "i:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) bam=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	r ) reference=$OPTARG ;;
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

	#Location of input FASTA files
        fasta=$directory/ribose-seq/results/$reference/$samples/Nucleotide-Frequencies/Ribonucleotides/$samples.fasta

	#OUTPUT
	#Location of output "ribose-seq" alignment directory
	output=$directory/ribose-seq/results/$reference/$samples/Nucleotide-Frequencies/Ribonucleotides

	#Create directory for output if it does not already exist
	if [[ ! -d $output ]];
	then
    		mkdir -p $output
	fi

	#Location of output files
	bed=$output/$samples.bed
	sam=$output/$samples.sam
	coordinates=$output/$samples.coordinates.bed
	positionsPositive0=$output/$samples.rNMPs.positive.0-based.txt
	positionsNegative0=$output/$samples.rNMPs.negative.0-based.txt
	positionsPositive1=$output/$samples.rNMPs.positive.1-based.txt
	positionsNegative1=$output/$samples.rNMPs.negative.1-based.txt
	
	#COORDINATES (0-BASED) of SEQUENCING READS

	#Covert file from BAM to BED format using BEDtools software
	bedtools bamtobed -i $input > $bed

	#Convert file from BAM to SAM format using SAMtools software
	samtools view $input > $sam

	#Extract read coordinates, sequences, and strands from BED and SAM files and save it to new file
	paste $bed $sam $fasta | awk -v "OFS=\t" '{print $1, $2, $3, $4, $16, $6, $21}' > $coordinates

	#0-BASED COORDINATES OF rNMPs

	#Obtain positions of rNMPs (3’ end of each mapped read) for positive strand:
	bedtools genomecov -3 -strand + -bg -ibam $input > $positionsPositive0

	#Obtain positions of rNMPs (3’ end of each mapped read) for negative strand:
	bedtools genomecov -3 -strand - -bg -ibam $input > $positionsNegative0

	#1-BASED COORDINATES OF	rNMPs

	#Obtain positions of rNMPs (3’ end of each mapped read) for positive strand:
	bedtools genomecov -3 -strand + -d -ibam $input > $positionsPositive1

	#Remove rows where genome coverage equals 0
	awk '$3 != 0' $positionsPositive1 > temporary

	#Change filename back to original
	mv temporary $positionsPositive1

	#Obtain positions of rNMPs (3’ end of each mapped read) for negative strand:
	bedtools genomecov -3 -strand - -d -ibam $input > $positionsNegative1

	#Remove rows where genome coverage equals 0
	awk '$3 != 0' $positionsNegative1 > temporary

	#Change filename back to original
	mv temporary $positionsNegative1

done

mkdir $output/tables
