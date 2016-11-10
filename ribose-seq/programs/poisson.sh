#!/usr/bin/env bash

#Author: Alli Gombolay
#This program bins rNMPs into 2.5kb windows in reference genome

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 4_Poisson.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of genome (sacCer2, nuclear, chrM, etc.)
	-r Reference genome assembly version (sacCer2, etc.)
	-d Local directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Use getopts function to create the command-line options ([-i], [-s], [-r], [-d], and [-h])
while getopts "i:s:r:d:h" opt; do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Separate genome into 2.5 kb windows
bedtools makewindows -g chrom.sizes -w 2500 > genome.windows.bed

#Sort BED file of rNMP coordinates (sort first and second columns of file)
sort -k1,1 -k2,2n FS15.sorted.bed > temporary && mv temporary FS15.sorted.bed

#Determine intersect regions of the two BED files and then bin data
bedtools intersect -a genome.windows.bed -b FS15.sorted.bed -c -sorted
