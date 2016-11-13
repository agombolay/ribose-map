#!/usr/bin/env bash
#Author: Alli Gombolay
#This program determines the coordinates of the rNMPs (3' position of aligned reads)

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 2_rNMP-Coordinates.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of genome (sacCer2, nuclear, chrM, etc.)
	-r Reference genome assembly version (sacCer2, etc.)
	-d Local directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Use getopts function to create the command-line options ([-i], [-s], [-r], [-d], and [-h])
while getopts "i:s:r:d:h" opt; do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Determine rNMP coordinates for each sample
for sample in ${sample[@]}; do

#############################################################################################################################
	#Input/Output
	
	#Location of input file
	bam=$directory/ribose-seq/results/$reference/$sample/Alignment/$sample.bam
	
	#Location of output directory
	output=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset

	#Create directory if it does not already exist
	mkdir -p $output
	
	#Remove any older versions of the output files
	rm -f $output/*.txt
	
	#Location of output files	
	fastq=$output/$sample.aligned-reads.fq
	fasta=$output/$sample.aligned-reads.fa
	sequences=$output/$sample.sequences.txt
	
	bed=$output/$sample.aligned-reads.bed
	coverage=$output/$sample.rNMP-coverage.txt
	readInformation=$output/$sample.read-information.bed
	
	riboCoordinates=$output/$sample.rNMP-coordinates.bed

#############################################################################################################################
	#STEP 1: Covert BAM alignment file to FASTA format

	#Convert BAM file to FASTQ
	samtools bam2fq $bam > $fastq

	#Convert FASTQ file to FASTA
	seqtk seq -A $fastq > $fasta
	
	#Extract sequences from FASTA file
	#grep -v '>' $fasta > temporary && mv temporary $fasta
	grep -v '>' $fasta > $sequences

#############################################################################################################################
	#STEP 2: Obtain rNMP coordinates from aligned reads

	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > $bed
	
	#Obtain coverage of 3' positions of reads
	bedtools genomecov -3 -bg -ibam $bam > $coverage
	
	#Extract read coordinates, sequences, and strands from BED and FASTA files
	paste $bed $sequences | awk -v "OFS=\t" '{print $1, $2, $3, $4, $6, $7}' > $readInformation
	
	#Determine rNMP coordinates from reads aligned to positive strand of DNA
	positiveReads=$(awk -v "OFS=\t" '$5 == "+" {print $1, ($3 - 1), $3, " ", " ", $5}' $readInformation)

	#Determine rNMP coordinates from reads aligned to negative strand of DNA
	negativeReads=$(awk -v "OFS=\t" '$5 == "-" {print $1, $2, ($2 + 1), " ", " ", $5}' $readInformation)

	#Select only rNMP coordinates located in nuclear DNA
	if [ $subset == "nuclear" ]; then
		echo -e "\v$positiveReads \v$negativeReads" | grep -v 'chrM' - > $riboCoordinates
	#Select only rNMP coordinates located in mitochondrial DNA
	elif [ $subset == "chrM" ]; then
		echo -e "\v$positiveReads \v$negativeReads" | grep 'chrM' - > $riboCoordinates
	#Select all rNMP coordinates located in genomic DNA
	else
		echo -e "\v$positiveReads \v$negativeReads" > $riboCoordinates
	fi
	
done
