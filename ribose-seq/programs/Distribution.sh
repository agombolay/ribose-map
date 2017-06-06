#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program counts number of positions in nucleus and mitochondria with 0...X rNMPs

#Usage statement
function usage () {
	echo "Usage: Distribution.sh [-s] 'Sample(s)' [-r] 'Reference' [-d] 'Directory' [-h]
	-s Sample name(s) (e.g., FS1, FS2, FS3)
	-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
	-d Directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:r:d:h" opt; do
    case $opt in
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
	for subset in "mito"; do

#############################################################################################################################
	#Create directory
	mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Distribution
	
	#Input files
	bed=$directory/Ribose-Map/References/$reference.bed
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-Coordinates.$subset.bed
	#bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample-MappedReads.bam

	#Output file
	dataset=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-ObservedCounts.$subset.txt

	#Remove old files
	rm -f $dataset temp{1..3}.txt windows.bed
#############################################################################################################################
	
	#STEP 1: Divide genome into windows and count number of rNMPs in each window
	
	#Divide chromosomes of reference into windows
	if [ $subset == "all" ]; then
		bedtools makewindows -g $bed -w 1 > windows.bed
	elif [ $subset == "mito" ]; then
		bedtools makewindows -g $bed -w 1 | grep -E '(chrM|MT)' > windows.bed
	elif [ $subset == "nucleus" ]; then
		bedtools makewindows -g $bed -w 1 | grep -vE '(chrM|MT)' > windows.bed
	fi

	#Determine regions of BED files that intersect and count number of intersections
	bedtools intersect -a windows.bed -b $coordinates -c -sorted -nonamecheck > temp1.txt
	#bedtools genomecov -d -3 -ibam $bam > temp1.txt
	
	#for subset in "all" "mito" "nucleus"; do
	
	#Output file
	#dataset=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-ObservedCounts.$subset.txt

	#Remove old file
	#rm -f $directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-ObservedCounts.$subset.txt
	
	#Divide chromosomes of reference into windows
	#if [ $subset == "all" ]; then
	#	cat temp1.txt > temp2.txt
	#elif [ $subset == "mito" ]; then
	#	grep -E '(chrM|MT)' temp1.txt > temp2.txt
	#elif [ $subset == "nucleus" ]; then
	#	grep -vE '(chrM|MT)' temp1.txt > temp2.txt
	#fi
	
	#Sort by # of rNMPs
	sort -k4n temp1.txt > temp2.txt

	#Maximum # of rNMPs in observed data
	max=$(tail -1 temp2.txt | awk '{print $4}' -)

	#Number of positions containing 0...max # of rNMPs
	for i in $(seq 0 $max); do
		awk '$4 == ('$i')' temp2.txt | wc -l >> temp3.txt
	done

#############################################################################################################################
	#STEP 6: Create and save dataset file containing observed counts of rNMPs
	
	#Add column names to header line
	echo -e "rNMPs\tPositions" > $dataset

	#Add number of positions containing 0...max # of rNMPs
	paste <(echo "$(seq 0 $max)") <(cat temp3.txt) >> $dataset

	#Print completion status
	echo "Observed counts for $sample ($subset) have been determined"
	
	#Remove temp files
	#rm -f temp{1..3}.txt windows.bed

	done
done
#rm -f temp1.txt
