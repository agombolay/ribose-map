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

#Determine coordinates
for sample in ${sample[@]}; do

#############################################################################################################################
		#Create directory
		mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Hotspots
		
		#Input file
		coverage=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Coverage.$subset.bed

#############################################################################################################################
		for subset in "mito" "nucleus"; do
			for i in "I" "II"
			sed 's/chrI/chr1/' $coverage
			sed 's/chrII/chr2/' $coverage
			sed 's/chrIII/chr3/' $coverage
			sed 's/chrIV/chr4/' $coverage
			sed 's/chrV/chr5/' $coverage
			sed 's/chrVI/chr6/' $coverage
			sed 's/chrVII/chr7/' $coverage
			sed 's/chrVIII/chr8/' $coverage
			sed 's/chrIX/chr9/' $coverage
			sed 's/chrX/chr10/' $coverage
			sed 's/chrXI/chr11/' $coverage
			sed 's/chrXII/chr12/' $coverage
			sed 's/chrXIII/chr13/' $coverage
			sed 's/chrXIV/chr14/' $coverage
			sed 's/chrXV/chr15/' $coverage
			sed 's/chrXVI/chr16/' $coverage
			sed 's/chrXVII/chr17/' $coverage
			sed 's/chrXVIII/chr18/' $coverage
			sed 's/chrXIX/chr19/' $coverage
			sed 's/chrXX/chr20/' $coverage
			sed 's/chrXXI/chr21/' $coverage
			sed 's/chrXXII/chr22/' $coverage
		done

done
