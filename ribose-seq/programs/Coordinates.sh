#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Determine the genomic coordinates of rNMPs
#2. Can be applied to rNMP sequencing technique

#Usage statement
function usage () {
	echo "Usage: Coordinates.sh [options]
	-d Ribose-Map repository
	-s Name of sequenced library
	-i Basename of Bowtie2 index
	-t rNMP sequencing technique"
}

#Command-line options
while getopts "s:t:r:d:h" opt; do
    case $opt in
        s ) sample=$OPTARG ;;
	t ) technique=$OPTARG ;;
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

#############################################################################################################################
#Output directory
output=$directory/results/$reference/$sample/coordinates

#Input alignment file
bam=$directory/results/$reference/$sample/alignment/$sample.bam
	
#Create directory and remove old files
mkdir -p $output; rm -f $output/*.{bed}
		
#############################################################################################################################
#Determine coordinates for each technique
if [[ "$technique" == "ribose-seq" ]]; then
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,($3 - 1),$3," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,$2,($2 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
elif [[ "$technique" == "emRiboSeq" ]]; then
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$4 == "-" {print $1,$3,($3 + 1)," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$4 == "+" {print $1,($2 - 1),$2," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
elif [[ "$technique" == "HydEn-seq" ]] || [[ "Pu-seq" ]]; then
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$4 == "+" {print $1,($2 - 1),$2," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$4 == "-" {print $1,$3,($3 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
fi
	
#Sort chromosome coordinates of rNMPs
sort -k1,1V -k2,2n $output/temp2.bed > $output/$sample-Coordinates.bed

#############################################################################################################################
#Print completion status
echo "Status: Genomic coordinates of rNMPs were identified for $sample"
	
#Remove temporary files
rm -f $output/temp{1..2}.bed
