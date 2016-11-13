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
	rm -f $output/{*.txt,*.bed,*.fa,*.fq}
	
	#Location of output files	
	sequences=$output/$sample.sequences.txt; bed=$output/$sample.aligned-reads.bed;
	reads=$output/$sample.read-information.bed; coordinates=$output/$sample.rNMP-coordinates.bed

#############################################################################################################################
	#STEP 1: Extract sequences from BAM alignment file

	#Convert BAM to FASTA file then extract sequences from FASTA
	samtools bam2fq $bam | seqtk seq -A - | grep -v '>' - > $sequences

#############################################################################################################################
	#STEP 2: Obtain rNMP coordinates from aligned reads

	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > $bed
	
	#Extract read coordinates, sequences, and strands from BED and FASTA files
	paste $bed $sequences | awk -v "OFS=\t" '{print $1, $2, $3, $4, $6, $7}' > $reads
	
	#Determine rNMP coordinates from reads aligned to positive strand of DNA
	#awk -v "OFS=\t" '$5 == "+" {print $1, ($3 - 1), $3, " ", " ", $5}' $reads > positiveReads

	#Determine rNMP coordinates from reads aligned to negative strand of DNA
	#awk -v "OFS=\t" '$5 == "-" {print $1, $2, ($2 + 1), " ", " ", $5}' $reads > negativeReads

	#Obtain rNMP coordinates
	if [ $subset == "nuclear" ]; then
		#Select only rNMP coordinates located in nuclear DNA
		echo "$positiveReads $negativeReads" | grep -v 'chrM' - > $coordinates
	elif [ $subset == "chrM" ]; then
		#Select only rNMP coordinates located in mitochondrial DNA
		echo "$positiveReads $negativeReads" | grep 'chrM' - > $coordinates
	else
		#Select all rNMP coordinates located in genomic DNA
		#echo "$positiveReads $negativeReads" > $coordinates
		#cat positiveReads negativeReads > $coordinates
		positiveReads=$(awk -v "OFS=\t" '$5 == "+" {print $1, ($3 - 1), $3, " ", " ", $5}' $reads)
		negativeReads=$(awk -v "OFS=\t" '$5 == "-" {print $1, $2, ($2 + 1), " ", " ", $5}' $reads)
		awk -v "OFS=\t" '{printf $positiveReads}' > file.bed
	fi
	
done
