#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program determines the coordinates of rNMPs (3' position of aligned reads)

#Usage statement
function usage () {
	echo "Usage: Coordinates.sh [-i] 'Sample(s)' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Sample name(s) (FS1, FS2, FS3 etc.)
	-r Reference genome (sacCer2, pombe, ecoli, mm9, hg38, etc.)
	-d Directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Command-line options
while getopts "i:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Allow only one input argument
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

subset=("genome" "nucleus" "mitochondria")

#Determine coordinates
for sample in ${sample[@]}; do
	for subset in ${subset[@]}; do

#############################################################################################################################
		#Input file
		bam=$directory/ribose-seq/results/$reference/$sample/Alignment/$sample.bam
	
		#Output directory
		output=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset

		#Create directory
		mkdir -p $output
	
		#Remove older versions of files
		rm -f $output/{*.txt,*.bed,*.fa,*.fq}
	
		#Output files
		reads=$output/$sample.read-information.$subset.txt
		coordinates=$output/$sample.rNMP-coordinates.$subset.bed

#############################################################################################################################
		#STEP 1: Extract sequences from BAM alignment file

		#Convert BAM to FASTA file then extract sequences from FASTA
		samtools bam2fq $bam | seqtk seq -A - | grep -v '>' - > temp1

#############################################################################################################################
		#STEP 2: Obtain rNMP coordinates from aligned reads

		#Covert BAM file to BED format
		bedtools bamtobed -i $bam > temp2
		#Extract read coordinates, sequences, and strand information
		paste temp2 temp1 | awk -v "OFS=\t" '{print $1, $2, $3, $4, $6, $7}' > $reads
		#Obtain coordinates of rNMPs located on positive strand of DNA
		positiveReads=$(awk -v "OFS=\t" '$5 == "+" {print $1, ($3 - 1), $3, " ", " ", $5}' $reads)
		#Obtain coordinates of rNMPs located on negative strand of DNA
		negativeReads=$(awk -v "OFS=\t" '$5 == "-" {print $1, $2, ($2 + 1), " ", " ", $5}' $reads)
	
		#Combine +/- coordinates into one file for processing
		cat <(echo "$positiveReads") <(echo "$negativeReads") > temp3
		
		#Subset and/or sort coordinates
		if [ $subset == "genome" ]; then
			sort -k1,1 -k2,2n temp3 > $coordinates
		elif [ $subset == "mitochondria" ]; then
			grep -E '(chrM|MT)' temp3 | sort -k1,1 -k2,2n - > $coordinates
		elif [ $subset == "nucleus" ]; then
			grep -v -E '(chrM|MT*|AB*|chrEBV|chrUN*|*random)' temp3 | sort -k1,1 -k2,2n - > $coordinates
		fi
	done
done

#Remove temp files
rm temp1 temp2 temp3
