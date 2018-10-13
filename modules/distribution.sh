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
output=$repository/results/$sample/distribution-$quality
rm -rf $output; mkdir -p $output

#############################################################################################################################
#Create FASTA index file and BED file for reference
samtools faidx $fasta && cut -f 1,2 $fasta.fai > $output/reference.bed

if [[ ! $technique ]]; then
	#Create file of normalized rNMP counts if user did not use Coordinate Module
	cut -f1,2,3,6 "$2" | uniq -c - | mawk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' > $output/$sample.counts.tab
	mawk -v "OFS=\t" -v total="$(wc -l < "$2")" '{print $1, $2, $3, $4, $5/total*100}' $output/$sample.counts.tab > $output/$sample.normalized.tab
fi

#Save coverage of rNMPs per chromosome to separate files
for chromosome in $( awk '{print $1}' $output/reference.bed ); do
	
	if [[ $(grep -w "$chromosome" $repository/results/$sample/coordinate-$quality/$sample.normalized.tab | wc -l) > 0 ]]; then
		grep -w "$chromosome" $repository/results/$sample/coordinate-$quality/$sample.normalized.tab > $output/$sample-$chromosome.tab
	fi
	
done

if [[ $(wc -l < $repository/results/$sample/coordinate-$quality/$sample.counts.tab) > 0 ]]; then

	#Add trackline for forward strand to input into UCSC genome browser
	echo "track type=bedGraph name="$sample-ForwardStrand" description="$sample-ForwardStrand" color=0,128,0 visibility=full" > $output/$sample-Forward.bg
		
	#Add trackline for reverse strand to input into UCSC genome browser
	echo "track type=bedGraph name="$sample-ReverseStrand" description="$sample-ReverseStrand" color=0,0,255 visibility=full" > $output/$sample-Reverse.bg
		
	#Rearrange forward strand file so format is the same as bedgraph format
	awk -v "OFS=\t" '$4 == "+" {print $1, $2, $3, $5}' $repository/results/$sample/coordinate-$quality/$sample.counts.tab >> $output/$sample-Forward.bg

	#Rearrange reverse strand file so format is the same as bedgraph format
	awk -v "OFS=\t" '$4 == "-" {print $1, $2, $3, $5}' $repository/results/$sample/coordinate-$quality/$sample.counts.tab >> $output/$sample-Reverse.bg

fi
#############################################################################################################################
#Print status
echo "Status: Distribution Module for $sample is complete"
