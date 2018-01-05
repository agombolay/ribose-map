#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program determines the coordinates of rNMPs
#Ribose-seq rNMP = RC of 5' position of aligned reads

#Usage statement
function usage () {
	echo "Usage: Coordinates.sh [options]
	-s Sample name(s) (e.g., FS1, FS2, FS3)
	-r Reference genome/Basename of Bowtie2 index
	-d Ribose-Map directory (e.g., /path/to/Ribose-Map)"
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
	bam=$directory/Results/$reference/$sample/Alignment/$sample.bam

	#Output directory
	output=$directory/Results/$reference/$sample/Coordinates

#############################################################################################################################
	#Create directory
	mkdir -p $output

	#Remove any old files
	rm -f $output/*.{txt,bed}
		
#############################################################################################################################
	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > $output/temp1.txt
		
	#Extract read coordinates, sequences, and strand information
	awk -v "OFS=\t" '{print $1, $2, $3, $6}' $output/temp1.txt > $output/temp2.txt
		
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$4 == "-" {print $1,($3 - 1),$3," "," ",$4}' $output/temp2.txt > $output/temp3.txt 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$4 == "+" {print $1,$2,($2 + 1)," "," ",$4}' $output/temp2.txt > $output/temp4.txt

	#Combine and save +/- coordinates into one file for later
	cat $output/temp3.txt $output/temp4.txt | sort -k1,1V -k2,2n > $output/$sample-Coordinates.bed

#############################################################################################################################
	#Remove temp files
	rm -f $output/temp{1..4}.txt

	#Print completion status
	echo "Coordinates of rNMPs for $sample have been determined"
	
done
