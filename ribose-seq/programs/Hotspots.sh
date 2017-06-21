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
	for subset in "mito" "nucleus"; do
	
		#Create directory
		mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Hotspots
		
		#Input file
		coverage=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Coverage.$subset.bed

		#Output file
		hotspots=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Hotspots.$subset.bed
		
		sed 's/chrI/chr1/' $coverage | sort -k1 > $hotspots
		sed 's/chrII/chr2/' $coverage | sort -k1 > $hotspots
		sed 's/chrIII/chr3/' $coverage | sort -k1 > $hotspots
		sed 's/chrIV/chr4/' $coverage | sort -k1 > $hotspots
		sed 's/chrV/chr5/' $coverage | sort -k1 > $hotspots
		sed 's/chrVI/chr6/' $coverage | sort -k1 > $hotspots
		sed 's/chrVII/chr7/' $coverage | sort -k1 > $hotspots
		sed 's/chrVIII/chr8/' $coverage | sort -k1 > $hotspots
		sed 's/chrIX/chr9/' $coverage | sort -k1 > $hotspots
		sed 's/chrX/chr10/' $coverage | sort -k1 > $hotspots
		sed 's/chrXI/chr11/' $coverage | sort -k1 > $hotspots
		sed 's/chrXII/chr12/' $coverage | sort -k1 > $hotspots
		sed 's/chrXIII/chr13/' $coverage | sort -k1 > $hotspots
		sed 's/chrXIV/chr14/' $coverage | sort -k1 > $hotspots
		sed 's/chrXV/chr15/' $coverage | sort -k1 > $hotspots
		sed 's/chrXVI/chr16/' $coverage | sort -k1 > $hotspots
		sed 's/chrXVII/chr17/' $coverage | sort -k1 > $hotspots
		sed 's/chrXVIII/chr18/' $coverage | sort -k1 > $hotspots
		sed 's/chrXIX/chr19/' $coverage | sort -k1 > $hotspots
		sed 's/chrXX/chr20/' $coverage | sort -k1 > $hotspots
		sed 's/chrXXI/chr21/' $coverage | sort -k1 > $hotspots
		sed 's/chrXXII/chr22/' $coverage | sort -k1 > $hotspots
	done
done
