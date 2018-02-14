#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Determine the chromosome coordinates of rNMPs
#2. Can be applied to any rNMP sequencing technique

. /data2/users/agombolay3/Ribose-Map/config.txt

output=$directory/results/$sample/coordinates
mkdir -p $output
			
#############################################################################################################################
#Determine coordinates for each technique
if [[ "$technique" == "ribose-seq" ]]; then
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $directory/results/$sample/alignment/$sample.bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,($3 - 1),$3," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,$2,($2 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
elif [[ "$technique" == "emRiboSeq" ]]; then
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,$3,($3 + 1)," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,($2 - 1),$2," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
elif [[ "$technique" == "HydEn-seq" ]] || [[ "Pu-seq" ]]; then
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,($2 - 1),$2," "," ","+"}' $output/temp1.bed > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,$3,($3 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp2.bed
	
fi
	
#Sort chromosome coordinates of rNMPs
sort -k1,1V -k2,2n $output/temp2.bed > $output/$sample-Coordinates.bed

#############################################################################################################################
#Print completion status
echo "Status: Genomic coordinates of rNMPs were identified for $sample"
	
#Remove temporary files
rm -f $output/temp{1..2}.bed
