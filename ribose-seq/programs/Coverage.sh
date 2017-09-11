#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates the coveage at each rNMP position

#Usage statement
function usage () {
	echo "Usage: Coverage.sh [options]
		-s Sample name(s) (e.g., FS1, FS2, FS3)
		-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository/Ribose-Map)"
}

#Command-line options
while getopts "s:r:d:h" opt; do
    case "$opt" in
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

#Determine coordinates
for sample in ${sample[@]}; do

#############################################################################################################################
	#Input file
	bed=$directory/References/$reference.bed
	bam=$directory/Results/$reference/$sample/Alignment/$sample.bam
	bedgraph1=$directory/Results/$reference/$sample/Hotspots/$sample-Forward.bedgraph
	bedgraph2=$directory/Results/$reference/$sample/Hotspots/$sample-Reverse.bedgraph
	
	#Output directory
	output=$directory/Results/$reference/$sample/Coverage
	
	#Create directory
	mkdir -p $output
#############################################################################################################################
	
	if [[ -s $bedgraph* ]]; then
		
		#Remove old files
		rm -f $output/$sample-*.bed

		#Combine forward and reverse bedgraph files
		cat $bedgraph1 bedgraph2 > $output/$sample-rNMPs.bed
		
		#Calculate coverage of each position in genome
		bedtools genomecov -dz -ibam $bam -g $bed > $output/$sample-genome.bed
		
		#Determine coverage at each rNMP position in genome
		bedtools intersect -a $output/$sample-rNMPs.bed -b $output/$sample.bed > $output/$sample.bed

	fi

#Remove temporary files
rm -f $output/$sample-rNMPs.bed $output/$sample-genome.bed