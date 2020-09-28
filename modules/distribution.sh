#!/usr/bin/env bash

#Author: Alli Lauren Gombolay
#Creates bedgraph files of per nucleotide rNMP coverage for each DNA strand

#############################################################################################################################

#Load config file
. "$1"

#Create output directory and remove any old files
output=$repository/results/$sample/distribution$quality
rm -r $output; mkdir -p $output

#############################################################################################################################

for region in $unit "chromosomes"; do

	#Calculate normalized counts of rNMPs
	mawk -v "OFS=\t" -v total="$(wc -l < $repository/results/$sample/coordinate$quality/$sample-$region.coords.bed)" '{print $1, $2, $3, $6, $7/total*100}' \
	$repository/results/$sample/coordinate$quality/$sample-$region.counts.tab > $repository/results/$sample/coordinate$quality/$sample-$region.normalized.tab

	#Save coverage of rNMPs per chromosome to separate files
	for unit in $( awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes ); do
	
		if [[ $(grep -w "$unit" $repository/results/$sample/coordinate$quality/$sample-$region.normalized.tab | wc -l) > 0 ]]; then
			grep -w "$unit" $repository/results/$sample/coordinate$quality/$sample-$region.normalized.tab > $output/$sample-$unit.tab
		fi
	
	done

	if [[ $(wc -l < $repository/results/$sample/coordinate$quality/$sample.counts.tab) > 0 ]]; then

		#Add trackline for forward strand to input into UCSC genome browser
		echo "track type=bedGraph name="$sample-ForwardStrand" description="$sample-ForwardStrand" color=0,128,0 visibility=full" > $output/$sample-Forward.bg
		
		#Add trackline for reverse strand to input into UCSC genome browser
		echo "track type=bedGraph name="$sample-ReverseStrand" description="$sample-ReverseStrand" color=0,0,255 visibility=full" > $output/$sample-Reverse.bg
		
		#Rearrange forward strand file so format is the same as bedgraph format
		awk -v "OFS=\t" '$4 == "+" {print $1, $2, $3, $5}' $repository/results/$sample/coordinate$quality/$sample.counts.tab >> $output/$sample-Forward.bg

		#Rearrange reverse strand file so format is the same as bedgraph format
		awk -v "OFS=\t" '$4 == "-" {print $1, $2, $3, $5}' $repository/results/$sample/coordinate$quality/$sample.counts.tab >> $output/$sample-Reverse.bg

	fi

done

#############################################################################################################################

#Print status
echo "Status: Distribution Module for $sample is complete"
