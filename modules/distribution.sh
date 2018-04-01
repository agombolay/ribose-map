#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Creates bedgraph files for forward and reverse strands
#2. Saves coverage of rNMPs per chromosome to separate files

#############################################################################################################################
#Load config file
. "$1"

#Create output directory and remove any old files
output=$repository/results/$sample/distribution; rm -rf $output; mkdir -p $output

#############################################################################################################################
#Create FASTA index and BED file for reference genome
samtools faidx $fasta && cut -f 1,2 $fasta.fai > $output/reference.bed

#Create file of rNMP coverage at chromosome coordinates
uniq -c $repository/results/$sample/coordinates/$sample.bed | awk -v "OFS=\t" '{print $2,$3,$4,$1,$5}' > $output/temp.tab

#############################################################################################################################
#Add trackline for forward strand to input into UCSC genome browser
echo "track type=bedGraph name="$sample-ForwardStrand" description="$sample-ForwardStrand" color=0,128,0 visibility=full" > $output/$sample-Forward.bg
		
#Add trackline for reverse strand to input into UCSC genome browser
echo "track type=bedGraph name="$sample-ReverseStrand" description="$sample-ReverseStrand" color=0,0,255 visibility=full" > $output/$sample-Reverse.bg
		
#Rearrange forward strand file so format is the same as bedgraph format
awk -v "OFS=\t" '$5 == "+" {print $1,$2,$3,$4}' $output/temp.tab >> $output/$sample-Forward.bg

#Rearrange reverse strand file so format is the same as bedgraph format
awk -v "OFS=\t" '$5 == "-" {print $1,$2,$3,$4}' $output/temp.tab >> $output/$sample-Reverse.bg

#############################################################################################################################
#Calculate normalized per-nucleotide coverage
for coverage in $(ls $output/temp.tab); do
	total=$(samtools view -c $repository/results/$sample/alignment/$sample.bam)
	awk -v "OFS=\t" -v total="$total" '{print $1,$2,$3,$4/total*1000000,$5}' $coverage > $output/normalized.bed
done

#Save coverage of rNMPs per chromosome to separate files
for chromosome in $( awk '{print $1}' $output/reference.bed ); do
	grep -w "$chromosome" $output/normalized.bed > $output/$sample-$chromosome.bed
done

#############################################################################################################################
#Print status
echo "Status: Distribution module for $sample is complete"
	
#Remove temporary files
rm $output/reference.bed $output/temp.tab $output/normalized.bed
