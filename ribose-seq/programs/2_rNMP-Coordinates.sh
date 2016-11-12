#!/usr/bin/env bash
#Author: Alli Gombolay
#This program determines the coordinates of the rNMPs (3' position of aligned reads)

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 5_Ribonucleotide-Frequencies.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
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

#Calculate nucleotide frequencies for each sample
for sample in ${sample[@]}; do

	#Location of "Alignment" directory
	directory1=$directory/ribose-seq/results/$reference/$sample/Alignment

	#Location of "Nucleotide-Frequencies" directory
	directory2=$directory/ribose-seq/results/$reference/$sample/rNMP-Coordinates

##########################################################################################################################################
	#STEP 1: Covert BAM alignment file to FASTA format

	#Location of input file
	bam=$directory1/$sample.bam

	#Location of output directory
	output=$directory2/$subset

	#Create directory if it does not already exist
	mkdir -p $output

	#Remove any older versions of the output files
	rm -f $output/*.txt
	
	#Location of output files	
	fastq=$output/$sample.aligned-reads.fq
	fasta=$output/$sample.aligned-reads.fa

	#Convert BAM file to FASTQ
	samtools bam2fq $bam > $fastq

	#Convert FASTQ file to FASTA
	seqtk seq -A $fastq > $fasta
	
	#Extract sequences from FASTA file
	grep -v '>' $fasta > temporary && mv temporary $fasta
	
##########################################################################################################################################
	#STEP 2: Obtain rNMP coordinates from aligned reads

	#Location of output files
	bed=$output/$sample.aligned-reads.bed
	coverage=$output/$sample.rNMP-coverage.0-based.txt
	
	readCoordinates=$output1/$sample.read-coordinates.bed
	readInformation=$output1/$sample.read-information.bed
	
	riboCoordinates1=$output1/$sample.rNMP-coordinates.0-based.genome.bed
	riboCoordinates2=$output1/$sample.rNMP-coordinates.0-based.$subset.bed
	
	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > $bed
	
	#Obtain genome coverage of aligned reads
	bedtools genomecov -3 -bg -ibam $bam > $coverage
	
	#Extract read coordinates, sequences, and strands from BED and SAM files
	paste $bed $fasta | awk -v "OFS=\t" '{print $1, $2, $3, $4, $6, $7}' > $readInformation
	
	#Determine rNMP coordinates from reads aligned to positive strand
	awk -v "OFS=\t" '$5 == "+" {print $1, ($3 - 1), $3, " ", " ", $5}' $readInformation > positive-reads.txt
	
	#Determine rNMP coordinates from reads aligned to negative strand
	awk -v "OFS=\t" '$5 == "-" {print $1, $2, ($2 + 1), " ", " ", $5}' $readInformation > negative-reads.txt
	
	#Combine reads on both strands (+ and -) into one file
	cat positive-reads.txt negative-reads.txt > $riboCoordinates1
	
	#Select only rNMP coordinates located in nuclear DNA
	if [ $subset == "nuclear" ]; then
		grep -v 'chrM' $riboCoordinates1 > $riboCoordinates2
	#Select only rNMP coordinates located in mitochondrial DNA
	elif [ $subset == "chrM" ]; then
		grep 'chrM' $riboCoordinates1 > $riboCoordinates2
	#Select all rNMP coordinates located in genomic DNA
	else
		cat $riboCoordinates1 > $riboCoordinates2
	fi
	
	#Remove intermediate files
	rm positive-reads.txt negative-reads.txt
	
done
