#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Determine the chromosome coordinates of rNMPs
#2. Can be applied to any rNMP sequencing technique

#############################################################################################################################
#Load config file
. "$1"

#Create output directory and remove any old files
output=$directory/results/$sample/coordinates; rm -rf $output; mkdir -p $output
			
#############################################################################################################################
if [[ read2 ]]; then
	samtools view -f67 -F260 $directory/results/$sample/alignment/$sample.bam | samtools sort - -o $output/temp.bam
	samtools index $output/temp.bam
fi

#Determine coordinates for each technique
if [[ "$technique" == "ribose-seq" ]]; then
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $directory/results/$sample/alignment/$sample.bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,($3 - 1),$3," "," ","+"}' $output/temp1.bed > $output/temp3.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,$2,($2 + 1)," "," ","-"}' $output/temp1.bed >> $output/temp3.bed
	
elif [[ "$technique" == "emRiboSeq" ]]; then
	
	#Create BED file for reference
	samtools faidx $reference && cut -f 1,2 $reference.fai > $output/reference.bed
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $directory/results/$sample/alignment/$sample.bam > $output/temp1.bed

	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,$3,($3 + 1)," "," ","+"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,($2 - 1),$2," "," ","-"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' >> $output/temp2.bed
	
	#Remove coordinates of rNMPs if the end position is greater than length of chromosome
	join -t $'\t' <(sort $output/reference.bed) <(sort $output/temp2.bed) | awk -v "OFS=\t" '$2 >= $4 { print $1,$3,$4," "," ",$5 }' > $output/temp3.bed
	
elif [[ "$technique" == "HydEn-seq" ]] || [[ "Pu-seq" ]]; then
	
	#Create fasta index and BED file for reference genome
	samtools faidx $reference && cut -f 1,2 $reference.fai > $output/reference.bed
	
	#Convert BAM file to BED format
	bedtools bamtobed -i $directory/results/$sample/alignment/$sample.bam > $output/temp1.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,($2 - 1),$2," "," ","+"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,$3,($3 + 1)," "," ","-"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' >> $output/temp2.bed

	#Remove coordinates of rNMPs if the end position is greater than length of chromosome
	join -t $'\t' <(sort $output/reference.bed) <(sort $output/temp2.bed) | awk -v "OFS=\t" '$2 >= $4 { print $1,$3,$4," "," ",$5 }' > $output/temp3.bed
fi
	
#Sort chromosome coordinates of rNMPs
sort -k1,1V -k2,2n $output/temp3.bed > $output/$sample.bed

#############################################################################################################################
#Print status
echo "Status: Coordinates module for $sample is complete"
	
#Remove temporary files
rm -f $output/reference.bed $output/temp{1..3}.bed
