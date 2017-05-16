#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program counts number of 2.5 kb windows in nucleus or mitochondria with 0...X rNMPs

#Usage statement
function usage () {
	echo "Usage: Distribution.sh [-i] 'Sample(s)' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Input sample(s) (e.g., FS1, FS2, FS3)
	-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
	-d Directory (e.g., /projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Command-line options
while getopts "i:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
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

#############################################################################################################################
#Input files
referenceBed=$directory/ribose-seq/reference/$reference.bed
coordinates=$directory/ribose-seq/results/$reference/$sample/Coordinates/genome/$sample-Coordinates.genome.bed

#Output directories
output1=$directory/ribose-seq/reference; output2=$directory/ribose-seq/results/$reference/$sample/Distribution/

#Create directories
mkdir -p $output1 $output2

#Output files
counts=$output2/$sample.observed-counts.txt
referenceWindows=$output1/$reference.windows.bed

#Separate chromosomes of reference into 2.5 kb windows
bedtools makewindows -g $referenceBed -w 2500 > $referenceWindows

#Determine regions of BED files that intersect and count number of intersections
#Remove rows where window size is < 2.5 kb and sort based on # of rNMPs in windows
bedtools intersect -a $referenceWindows -b $coordinates -c -sorted -nonamecheck \
| awk '{ $5 = $3 - $2 } 1' - | awk -v OFS='\t' '($5 == 2500 ) {print $1,$2,$3,$4}' - \
| sort -k4 -n - > temporary

#Maximum number of rNMPs in binned data
max=$(tail -1 temporary | awk '{print $4}' -)

#Determine number of windows with 0...maximum rNMPs
for i in $(seq 0 $max); do
	windows+=($(awk '$4 == ('$i')' temporary | wc -l))
done

#Add column names and # of windows with 0...maximum rNMPs
echo -e "rNMPs\tWindows" > $counts && paste <(echo "$(seq 0 $max)") \
<(cat <( IFS=$'\n'; echo "${windows[*]}" )) >> $counts

#Remove temporary file
rm temporary
