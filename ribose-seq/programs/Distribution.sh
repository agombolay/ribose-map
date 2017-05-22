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

subset=("all" "mito" "nucleus")

#Determine coordinates
for sample in ${sample[@]}; do
	for subset in "all" "mito" "nucleus"; do
#############################################################################################################################
#Input files
bed=$directory/Ribose-Map/Reference/$reference.bed
coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$subset/$sample-Coordinates.$subset.bed

#Output directories and files
output=$directory/Ribose-Map/Results/$reference/$sample/Distribution/; observed=$output/$sample-ObservedCounts.txt

#Create directory
mkdir -p $output

#Divide chromosomes of reference into windows
bedtools makewindows -g $bed -w $size > windows.bed

#Determine regions of BED files that intersect and count number of intersections
#Remove rows where window size is < size and sort based on # of rNMPs in windows
bedtools intersect -a windows.bed -b $coordinates -c -sorted -nonamecheck > temp1.txt

awk '{$5=$3-$2} 1' temp1.txt | awk -v OFS='\t' '($5=='$size') {print $1,$2,$3,$4}' | sort -k4n - > temp2.txt

#Maximum number of rNMPs in observed data
maximum=$(tail -1 temp2.txt | awk '{print $4}' -)
echo $maximum

#Determine # of windows with 0...maximum # of rNMPs
#for i in $(seq 0 $maximum); do
#	counts+=($(awk '$4 == ('$i')' temp2.txt | wc -l))
#done

#Add column names and # of windows with 0...maximum rNMPs
#echo -e "rNMPs\tWindows" > $counts && paste <(echo "$(seq 0 $max)") \
#<(cat <( IFS=$'\n'; echo "${windows[*]}" )) >> $observed

#Remove temp file
#rm temporary.txt

done
done
