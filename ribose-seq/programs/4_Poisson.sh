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

#Input directories
directory1=$directory/ribose-seq/reference/
directory2=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/rNMPs/$subset/

#Input files
referenceBED=$directory1/$reference.bed
riboCoordinates=$directory2/$sample.rNMP-coordinates.0-based.$subset.bed

#Output directories
output1=$directory/ribose-seq/reference/
output2=$directory/ribose-seq/results/$reference/$sample/Poisson

#Create directory if not present
mkdir -p $output1 $output2

#Output files
binnedData=$output2/$sample.binned.data.bed
referenceWindows=$output1/$reference.windows.bed
sortedBED=$output2/$sample.rNMP-coordinates.0-based.$subset.sorted.bed

#Separate reference genome into 2.5 kb (2,500 bp) windows
bedtools makewindows -g $referenceBED -w 2500 > $referenceWindows

#Sort BED file of ribonucleotide coordinates
sort -k1,1 -k2,2n $riboCoordinates > $sortedBED

#Determine regions of the two BED files that intersect with each other and then bin data
bedtools intersect -a $referenceWindows -b $sortedBED -c -sorted -nonamecheck > $binnedData

#Select only data of interest
if [ $subset == "nuclear" ]; then
	#For nuclear data, remove mitochondria data
	grep -v 'chrM' > temporary && mv temporary $binnedData
elif [ $subset == "chrM" ]; then
	#For nuclear data, remove nuclear data
	grep 'chrM' > temporary && mv temporary $binnedData
fi
