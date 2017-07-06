#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates the coveage at each rNMP position

#Usage statement
function usage () {
	echo "Usage: Hotspots.sh [options]
		-s Sample name(s) (e.g., FS1, FS2, FS3)
		-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
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

#Create directory
mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Hotspots

#Determine coordinates
for sample in ${sample[@]}; do

	bed=$directory/Ribose-Map/References/$reference.bed
	
	#Select chromosomes
	genome=$(awk '{print $1}' $bed)
	
	bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample.bam
	coverage=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Coverage.bed
	
	#Determine coverage at each 3' position
        bedtools genomecov -d -3 -ibam $bam > $coverage
	
	if [ -s $bam ]; then
	
        for chr in ${genome[@]}; do
		hotspots=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Hotspots.$chr.bed
        	grep -w "$chr" $coverage > $hotspots
	done
	
	#Print completion status
	echo "Hotspots in $sample have been located"
	
	fi
done
