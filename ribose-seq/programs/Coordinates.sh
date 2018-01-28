#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program determines the coordinates of rNMPs

#Usage statement
function usage () {
	echo "Usage: Coordinates.sh [options]
	-s Sample name(s) (e.g., FS1, FS2, FS3)
	-r Reference genome/Basename of Bowtie2 index
	-t rNMP sequencing technique used for library prep
	-d Ribose-Map directory (e.g., /path/to/Ribose-Map)"
}

#Command-line options
while getopts "s:t:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
        s ) sample=($OPTARG) ;;
	#Allow only one input argument
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

for sample in ${sample[@]}; do

#############################################################################################################################
	#Output directory
	output=$directory/Results/$reference/$sample/Coordinates

	#Path to input file
	bam=$directory/Results/$reference/$sample/Alignment/$sample.bam
	
#############################################################################################################################
	#Create directory
	mkdir -p $output

	#Remove old files
	rm -f $output/*.{bed}
		
#############################################################################################################################
	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > $output/temp1.bed
	
	#Determine coordinates for each sequencing technique
	if [[ "$technique" == "ribose-seq" ]]; then #ribose-seq
	
		#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
		awk -v "OFS=\t" '$6 == "-" {print $1,($3 - 1),$3," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
		#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
		awk -v "OFS=\t" '$6 == "+" {print $1,$2,($2 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
	elif [[ "$technique" == "emRiboSeq" ]]; then #emRiboSeq
	
		#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
		awk -v "OFS=\t" '$4 == "-" {print $1,$3,($3 + 1)," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
		#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
		awk -v "OFS=\t" '$4 == "+" {print $1,($2 - 1),$2," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
	elif [[ "$technique" == "HydEn-seq" ]]; then #HydEn-seq
	
		#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
		awk -v "OFS=\t" '$4 == "+" {print $1,($2 - 1),$2," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
		#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
		awk -v "OFS=\t" '$4 == "-" {print $1,$3,($3 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
	elif [[ "$technique" == "Pu-seq" ]]; then #Pu-seq
	
		#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
		awk -v "OFS=\t" '$4 == "+" {print $1,($2 - 1),$2," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
		#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
		awk -v "OFS=\t" '$4 == "-" {print $1,$3,($3 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
	fi
	
	#Sort chromosome coordinates of rNMPs
	sort -k1,1V -k2,2n $output/temp2.bed > $output/$sample-Coordinates.bed

#############################################################################################################################
	
	#Print completion status
	echo "Chromosome coordinates of rNMPs for $sample have been determined"
	
	#Remove temp files
	rm -f $output/temp{1..3}.bed
	
done
