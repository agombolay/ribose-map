#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates the coveage at each rNMP position

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [-s] 'Sample(s)' [-u] 'UMI' [-m] 'Min' [-p] 'Path' [-i] 'Index' [-d] 'Directory' [-h]
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
	for subset in "all" "mito" "nucleus"; do

#############################################################################################################################
		#Create directory
		mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Hotspots
		
		#Input file
		bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample-MappedReads.bam
		
		#Output file
		coverage=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Coverage.$subset.bed
		
		#Remove old files
		rm -f $directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Coverage.$subset.bed

#############################################################################################################################
		#Calculate coverage at each rNMP position
		bedtools genomecov -d -3 -ibam $bam > $coverage

done
done
