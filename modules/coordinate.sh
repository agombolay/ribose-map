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
output=$repository/results/$sample/coordinate-$quality; rm -rf $output; mkdir -p $output
			
#############################################################################################################################
#Convert alignment file to BED format
if [[ ! $read2 ]]; then
	#Remove unaligned reads
	samtools view -b -q $quality -@ $threads $repository/results/$sample/alignment/$sample.bam > $output/temp.bam

elif [[ $read2 ]]; then
	#Keep first read in pair
	samtools view -b -f67 -q $quality -@ $threads $repository/results/$sample/alignment/$sample.bam > $output/temp.bam
fi

#Convert BAM file to BED file
bedtools bamtobed -i $output/temp.bam > $output/temp1.bed

#Create FASTA index file and BED file for reference
samtools faidx $fasta && cut -f 1,2 $fasta.fai > $output/reference.bed
	
#Determine coordinates for each technique
if [[ $technique == "ribose-seq" ]]; then
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1, ($3 - 1), $3, $4, $5, "+"}' $output/temp1.bed > $output/temp2.bed
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1, $2, ($2 + 1), $4, $5, "-"}' $output/temp1.bed >> $output/temp2.bed
	
	#Sort coordinates by chromosome, position, and strand
	sort -k1,1 -k2,2n -k 6 $output/temp2.bed > $output/$sample.bed
	
elif [[ $technique == "emRiboSeq" ]]; then

	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1, $3, ($3 + 1), $4, $5, "+"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1, ($2 - 1), $2, $4, $5, "-"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' >> $output/temp2.bed

	#Remove coordinates of rNMPs if the end position is greater than length of chromosome
	#Sort by chromosome, position, and strand before joining files; otherwise, some data will be removed
	#Previous steps disorder temp2.bed since coordinates on forward/reverse strands are treated differently
	join -t $'\t' $output/reference.bed <(sort -k1,1 -k2,2n -k 6 $output/temp2.bed) | awk -v "OFS=\t" '$2 >= $4 { print $1, $3, $4, $5, $6, $7 }' > $output/$sample.bed
	
elif [[ $technique == "HydEn-seq" ]] || [[ $technique == "Pu-seq" ]]; then
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1, ($2 - 1), $2, $4, $5, "+"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1, $3, ($3 + 1), $4, $5, "-"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' >> $output/temp2.bed

	#Remove coordinates of rNMPs if the end position is greater than length of chromosome
	#Sort by chromosome, position, and strand before joining files; otherwise, some data will be removed
	#Previous steps disorder temp2.bed since coordinates on forward/reverse strands are treated differently
	join -t $'\t' $output/reference.bed <(sort -k1,1 -k2,2n -k 6 $output/temp2.bed) | awk -v "OFS=\t" '$2 >= $4 { print $1, $3, $4, $5, $6, $7 }' > $output/$sample.bed
fi

#Calculate raw and normalized (per 100) counts of rNMPs
total=$(wc -l < $output/$sample.bed)
cut -f1,2,3,6 $output/$sample.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' > $output/$sample.counts.tab
awk -v "OFS=\t" -v total="$total" '{print $1, $2, $3, $4, $5/total*100}' $output/$sample.counts.tab > $output/$sample.normalized.tab

#############################################################################################################################
#Remove temporary files
#rm -f $output/temp{1..2}.bed $output/temp.bam

#Print status
echo "Status: Coordinates Module for $sample is complete"
