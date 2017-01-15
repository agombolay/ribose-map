#!/usr/bin/env bash

#Author: Alli Gombolay
#This program counts the number of windows in the genome with 0...X number of rNMPs

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

#Input files
referenceBed=$directory/ribose-seq/reference/$reference.bed
sorted=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset/$sample.rNMP-coordinates.sorted.bed

#Output directories
output1=$directory/ribose-seq/reference/
output2=$directory/ribose-seq/results/$reference/$sample/Poisson

#Create directories
mkdir -p $output1 $output2

#Output files
binnedData=$output2/$sample.binned.data.bed
referenceWindows=$output1/$reference.windows.bed

counts1=$output2/$sample.probability-mass-function.counts.txt
counts2=$output2/$sample.cumulative-distribution.counts.txt

proportions1=$output2/$sample.probability-mass-function.proportions.txt
proportions2=$output2/$sample.cumulative-distribution.proportions.txt

#Separate reference genome into 2.5 kb windows
bedtools makewindows -g $referenceBed -w 2500 > $referenceWindows

#Select only data of interest
if [ $subset == "nuclear" ]; then
	#Determine regions of BED files that intersect and count number of overlaps (nuclear)
	bedtools intersect -a $referenceWindows -b $sorted -c -sorted -nonamecheck | grep -v 'chrM' -
	| awk '{ $5 = $3 - $2 } 1' - | awk '($5 == 2500 )' - > $binnedData
elif [ $subset == "chrM" ]; then
	#Determine regions of BED files that intersect and count number of overlaps (chrM)
	bedtools intersect -a $referenceWindows -b $sorted -c -sorted -nonamecheck | grep 'chrM' - > $binnedData
fi

#Remove rows where window size is < 2.5 kb
#awk '{ $5 = $3 - $2 } 1' $binnedData | awk '($5 == 2500 )' -

#Maximum value of genome coverage in BED file
maximum=$(sort -nk 4 $binnedData | tail -1 - | awk '{print $4}' -)

#Number of positions with X number of rNMPs
for i in $(seq 0 $maximum); do
	windows+=($(awk '$4 == ('$i')' $binnedData | wc -l))
done

#Probability Mass Function Counts
#Print number of windows in genomic region with exactly 0...X rNMPs
paste <(echo "$(seq 0 $maximum)") <(cat <( IFS=$'\n'; echo "${windows[*]}" )) > $counts1

#Cumulative Distribution Counts
#Print number of windows in genomic region with greater than or equal to 0...X rNMPs
for i in $(seq $(wc -l < $counts1) -1 1); do
	head -$(wc -l < $counts1) $counts1 | tail -${i} | awk '{ SUM += $2} END { print SUM }' >> $counts2
done

#Total number of windows
total=$(awk '{ SUM += $2} END { print SUM }' $counts1)

#Proportions of windows (P(X=x))
for i in ${windows[*]}; do
	values1+=($(echo "scale = 12; ($i/$total)" | bc | awk '{printf "%.12f\n", $0}'))
done

#Proportions of windows (P(X>=x))
for i in ${values1[*]}; do
	values2+=($(echo "scale = 12; (1-$i)" | bc | awk '{printf "%.12f\n", $0}'))
done

#Print proportions of windows (probability mass function and cumulative distribution)
paste <(echo "$(seq 0 $maximum)") <(cat <( IFS=$'\n'; echo "${values1[*]}" )) > $proportions1
paste <(echo "$(seq 0 $maximum)") <(cat <( IFS=$'\n'; echo "${values2[*]}" )) > $proportions2

#Calculate lambda value (mean) for Poisson Distribution
echo "Lambda:" $(echo "scale = 12; $(wc -l < $sorted)/$total" | bc | awk '{printf "%.5f\n", $0}')
