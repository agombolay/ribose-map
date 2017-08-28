#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program determines the coordinates of rNMPs (3' position of aligned reads)

#Usage statement
function usage () {
	echo "Usage: Coordinates.sh [options]
	-s Sample name(s) (e.g., FS1, FS2, FS3)
	-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
	-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
        s ) sample=($OPTARG) ;;
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

for sample in ${sample[@]}; do

#############################################################################################################################
	#Input file
	bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample.bam

	#Output directory
	output=$directory/Ribose-Map/Results/$reference/$sample/Coordinates
		
	#Create directory
	mkdir -p
		
#############################################################################################################################
	if [ -s $bam ]; then
		
		#Remove old file
		rm -f $output/$sample-ReadInfo.txt

		#Covert BAM file to BED format
		bedtools bamtobed -i $bam > $output/temp1.txt
		
		#Convert BAM file to FASTA then extract read sequences
		samtools bam2fq $bam | seqtk seq -A - | grep -v '>' - > $output/temp2.txt
		
		#Extract read coordinates, sequences, and strand information
		paste temp1.txt temp2.txt | awk -v "OFS=\t" '{print $1,$2,$3,$4,$6,$7}' > $output/$sample-ReadInfo.txt
		
		#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
		positiveReads=$(awk -v "OFS=\t" '$5 == "+" {print $1,($3 - 1),$3," "," ",$5}' $output/$sample-ReadInfo.txt)
			
		#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
		negativeReads=$(awk -v "OFS=\t" '$5 == "-" {print $1,$2,($2 + 1)," "," ",$5}' $output/$sample-ReadInfo.txt)
	
		#Combine and save +/- coordinates into one file for later
		cat <(echo "$positiveReads") <(echo "$negativeReads") > $output/temp3.txt

		#Sort coordinates
		#RomanNumerals.sh
		
		#sort -k1,1V -k2,2n $output/temp4.txt > $output/$sample-Coordinates.bed
		
	fi

#############################################################################################################################
	#Print completion status
	echo "Coordinates of rNMPs for $sample have been determined"
		
	#Remove temp files
	#rm -f $output/temp{1..3}.txt
		
done
