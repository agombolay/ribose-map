#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
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
output1=$directory/ribose-seq/reference; output2=$directory/ribose-seq/results/$reference/$sample/Poisson

#Create directories
mkdir -p $output1 $output2

#Output files
binnedData=$output2/$sample.binned.data.bed; referenceWindows=$output1/$reference.windows.bed
counts=$output2/$sample.PDF-windows.txt; proportions=$output2/$sample.PDF-proportions.txt

#Separate reference genome into 2.5 kb windows
bedtools makewindows -g $referenceBed -w 2500 > $referenceWindows

#Select only data of interest
if [ $subset == "nuclear" ]; then
	#Determine regions of BED files that intersect and count number of overlaps (nuclear)
	bedtools intersect -a $referenceWindows -b $sorted -c -sorted -nonamecheck | grep -v 'chrM' - | \
	#Remove rows where window size is < 2.5 kb
	awk '{ $5 = $3 - $2 } 1' - | awk '($5 == 2500 ) {print $1,$2,$3,$4}' - | sort -k4 -n - > $binnedData
elif [ $subset == "chrM" ]; then
	#Determine regions of BED files that intersect and count number of overlaps (chrM)
	bedtools intersect -a $referenceWindows -b $sorted -c -sorted -nonamecheck | grep 'chrM' - | \
	#Remove rows where window size is < 2.5 kb
	awk '{ $5 = $3 - $2 } 1' - | awk '($5 == 2500 ) {print $1,$2,$3,$4}' - | sort -k4 -n - > $binnedData
fi

#Maximum value of genome coverage in BED file
maximum=$(sort -nk 4 $binnedData | tail -1 - | awk '{print $4}' -)

#Determine number of windows with 0...maximum rNMPs
for i in $(seq 0 $maximum); do
	windows+=($(awk '$4 == ('$i')' $binnedData | wc -l))
done

#Print number of windows with 0...maximum rNMPs and save to file (input into R)
paste <(echo "$(seq 0 $maximum)") <(cat <( IFS=$'\n'; echo "${windows[*]}" )) > $counts
