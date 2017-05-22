#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program counts number of 2.5 kb windows in nucleus or mitochondria with 0...X rNMPs

#Usage statement
function usage () {
	echo "Usage: Distribution.sh [-s] 'Size' [-i] 'Sample(s)' [-r] 'Reference' [-d] 'Directory' [-h]
	-s Size of genomic windows
	-i Input sample(s) (e.g., FS1, FS2, FS3)
	-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
	-d Directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:i:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
        s ) size=($OPTARG) ;;
	i ) sample=($OPTARG) ;;
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
	for subset in "all"  "mito" "nucleus"; do
	
#############################################################################################################################
	#Input files
	bed=$directory/Ribose-Map/Reference/$reference.bed
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$subset/$sample-Coordinates.$subset.bed

	#Output directories and files
	output=$directory/Ribose-Map/Results/$reference/$sample/Distribution/; dataset=$output/$sample-Counts.$subset.txt

	#Create directory and remove old file
	mkdir -p $output; rm -f $output/$dataset
#############################################################################################################################
	
	#STEP 1: Divide genome into windows and count number of rNMPs in each window
	
	#Divide chromosomes of reference into windows
	if [ $subset == "all" ]; then
		bedtools makewindows -g $bed -w $size > windows.bed
	elif [ $subset == "mito" ]; then
		bedtools makewindows -g $bed -w $size | grep -E '(chrM|MT)' > windows.bed
	elif [ $subset == "nucleus" ]; then
		bedtools makewindows -g $bed -w $size | grep -vE '(chrM|MT)' > windows.bed
	fi

	#Determine regions of BED files that intersect and count number of intersections
	bedtools intersect -a windows.bed -b $coordinates -c -sorted -nonamecheck > temp1.txt

	#Remove rows where window size is greater than specified size and sort based on # of rNMPs in windows
	awk '{$5=$3-$2} 1' temp1.txt | awk -v OFS='\t' '($5=='$size') {print $1,$2,$3,$4}' | sort -k4n - > temp2.txt

	#Maximum # of rNMPs in observed data
	max=$(tail -1 temp2.txt | awk '{print $4}' -)

	#Number of windows containing 0...max # of rNMPs
	for i in $(seq 0 $max); do
		awk '$4 == ('$i')' temp2.txt | wc -l >> temp3.txt
	done

#############################################################################################################################
	#STEP 6: Create and save dataset file containing observed counts of rNMPs
	
	#Add column names to header line
	echo -e "rNMPs\tWindows" > $dataset

	#Add number of windows containing 0...max # of rNMPs
	paste <(echo "$(seq 0 $max)") <(cat temp3.txt) >> $dataset

	#Remove temp files
	rm temp{1..3}.txt

	done
done
