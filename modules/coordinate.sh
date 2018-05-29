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
output=$repository/results/$sample/coordinates; rm -rf $output; mkdir -p $output
			
#############################################################################################################################
#Convert alignment file to BED format
if [[ ! $read2 ]]; then
	#Remove unaligned reads
	samtools view -b -F2308 $repository/results/$sample/alignment/$sample.bam | bedtools bamtobed -i stdin > $output/temp1.bed

elif [[ $read2 ]]; then
	#Keep first read in pair
	samtools view -b -f67 -F2308 $repository/results/$sample/alignment/$sample.bam | bedtools bamtobed -i stdin > $output/temp1.bed
fi

#Determine coordinates for each technique
if [[ $technique == "ribose-seq" ]]; then
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,($3 - 1),$3,$4,$5,"+"}' $output/temp1.bed > $output/$sample.bed
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,$2,($2 + 1),$4,$5,"-"}' $output/temp1.bed >> $output/$sample.bed
	
elif [[ $technique == "emRiboSeq" ]]; then
	
	#Create FASTA index file and BED file for reference
	samtools faidx $fasta && cut -f 1,2 $fasta.fai > $output/reference.bed

	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,$3,($3 + 1),$4,$5,"+"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,($2 - 1),$2,$4,$5,"-"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' >> $output/temp2.bed

	#Remove coordinates of rNMPs if the end position is greater than length of chromosome
	join -t $'\t' $output/reference.bed $output/temp2.bed | awk -v "OFS=\t" '$2 >= $4 { print $1,$3,$4,$5,$6,$7 }' > $output/$sample.bed
	
elif [[ $technique == "HydEn-seq" ]] || [[ $technique == "Pu-seq" ]]; then
	
	#Create FASTA index file and BED file for reference
	samtools faidx $fasta && cut -f 1,2 $fasta.fai > $output/reference.bed
	
	#Obtain coordinates of rNMPs located on POSITIVE strand of DNA
	awk -v "OFS=\t" '$6 == "+" {print $1,($2 - 1),$2,$4,$5,"+"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' > $output/temp2.bed 
	
	#Obtain coordinates of rNMPs located on NEGATIVE strand of DNA
	awk -v "OFS=\t" '$6 == "-" {print $1,$3,($3 + 1),$4,$5,"-"}' $output/temp1.bed | awk -v "OFS=\t" '$2 >= 0 { print }' >> $output/temp2.bed

	#Remove coordinates of rNMPs if the end position is greater than length of chromosome
	join -t $'\t' $output/reference.bed $output/temp2.bed | awk -v "OFS=\t" '$2 >= $4 { print $1,$3,$4,$5,$6,$7 }' > $output/$sample.bed
fi

#Calculate per nucleotide coverage
total=$(wc -l < $output/$sample.bed)
cut -f1,2,3,6 $output/$sample.bed | uniq -c - | awk -v "OFS=\t" '{print $2,$3,$4,$5,$1}' > $output/$sample.counts.bed
awk -v "OFS=\t" -v total="$total" '{print $1,$2,$3,$4,$5/total*100}' $output/$sample.counts.bed > $output/$sample.normalized.bed

#############################################################################################################################
#Remove temporary files
#rm -f $output/reference.bed $output/temp{1..2}.bed

#Print status
echo "Status: Coordinates module for $sample is complete"
