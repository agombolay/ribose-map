#!/usr/bin/env bash

#Author: Alli Gombolay
#This program counts the number of rNMPs in each 2.5kb window of the reference genome

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

#############################################################################################################################

#Version 1: Proportion of positions that have x number of ribos

#Input files
bed=$directory/ribose-seq/reference/$reference.bed
bam=$directory/ribose-seq/results/$reference/$sample/Alignment/$sample.bam

#Output files
coverage=$directory/ribose-seq/results/$reference/$sample/Poisson/$sample.rNMP-coverage.bed
counts1=$directory/ribose-seq/results/$reference/$sample/Poisson/$sample.rNMP-counts.bed

#Obtain coverage at 3' positions of BAM file
if [ $subset == "nuclear" ]; then
	#Calculate total number of genome positions
	positions1=$(grep -v 'chrM' $bed | awk '{sum+=$2} END{print sum}' -)
	#Select only nuclear DNA regions from BED file
	bedtools genomecov -3 -bg -ibam $bam -g $bed | grep -v 'chrM' - > $coverage
elif [ $subset == "chrM" ]; then
	#Calculate total number of genome positions
	positions1=$(grep 'chrM' $bed | awk '{print $2}' -)
	#Select only mitochondrial DNA regions from BED file
	bedtools genomecov -3 -bg -ibam $bam -g $bed | grep 'chrM' - > $coverage
fi

#Maximum value of genome coverage in BED file
maximum=$(sort -nk 4 $coverage | tail -1 - | awk '{print $4}' -)

#Number of positions with X number of rNMPs
for i in $(seq 1 $maximum); do
	positions3+=($(awk '$4 == ('$i')' $coverage | wc -l))
done

#Number of positions with 0 rNMPs (total positions-positions with rNMPs)
positions2=$(echo "($positions1-$(wc -l $coverage | awk '{print $1}' -))" | bc)

#Print observed count data to fit to Poisson distribution
( IFS=$'\n'; echo -e "$positions2\n${positions3[*]}" ) > $counts1

#############################################################################################################################

#Version 2: Proportion of windows that have x number of ribos

#Input files
sorted=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset/$sample.rNMP-coordinates.sorted.bed

#Output directories
output1=$directory/ribose-seq/reference/
output2=$directory/ribose-seq/results/$reference/$sample/Poisson

#Create directory if not present
mkdir -p $output1 $output2

#Output files
windows=$output1/$reference.windows.bed
binned=$output2/$sample.binned.data.bed
counts2=$directory/ribose-seq/results/$reference/$sample/Poisson/$sample.rNMP-windows.bed

#Separate reference genome into 2.5 kb windows
bedtools makewindows -g $bed -w 2500 > $windows

#Obtain overlapping regions of BED files
if [ $subset == "nuclear" ]; then
	#Total number of windows in nuclear DNA
	windows1=$(grep -v 'chrM' $windows | wc -l)
	#Determine regions of BED files that intersect and count number of overlaps
	bedtools intersect -a $windows -b $sorted -c -sorted -nonamecheck | grep -v 'chrM' - > $binned
elif [ $subset == "chrM" ]; then
	#Total number of windows in mitochondrial DNA
	windows1=$(grep 'chrM' $windows | wc -l)
	#Determine regions of BED files that intersect and count number of overlaps
	bedtools intersect -a $windows -b $sorted -c -sorted -nonamecheck | grep 'chrM' - > $binned
fi

#Maximum value of rNMPs in a window of BED file
maximum=$(sort -nk 4 $binned | tail -1 - | awk '{print $4}' -)

#Number of windows with X number of rNMPs
for i in $(seq 1 $maximum); do
	counts1+=($(awk '$4 == ('$i')' $binned | wc -l))
done

#Number of windows with 0 rNMPs at any position in window
windows2=$(awk '$4 == 0' FS15.trimmed.v1.binned.data.bed | wc -l)

#Print observed count data to fit to Poisson distribution
( IFS=$'\n'; echo -e "$windows2\n${counts1[*]}" ) > $counts2
