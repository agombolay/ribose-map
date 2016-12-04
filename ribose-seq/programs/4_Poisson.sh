#!/usr/bin/env bash

#Author: Alli Gombolay
#This program counts the number of positions in the genome with 0...X number of rNMPs

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
bed=$directory/ribose-seq/reference/$reference.bed
bam=$directory/ribose-seq/results/$reference/$sample/Alignment/$sample.bam

#Output files
coverage=$directory/ribose-seq/results/$reference/$sample/Poisson/$sample.rNMP-coverage.bed
counts1=$directory/ribose-seq/results/$reference/$sample/Poisson/$sample.rNMP-counts.bed

#Obtain coverage at 3' positions of BAM file
if [ $subset == "nuclear" ]; then
	#Calculate total number of genome positions
	total=$(grep -v 'chrM' $bed | awk '{sum+=$2} END{print sum}' -)
	#Select only nuclear DNA regions from BED file
	bedtools genomecov -3 -bg -ibam $bam -g $bed | grep -v 'chrM' - > $coverage
elif [ $subset == "chrM" ]; then
	#Calculate total number of genome positions
	total=$(grep 'chrM' $bed | awk '{print $2}' -)
	#Select only mitochondrial DNA regions from BED file
	bedtools genomecov -3 -bg -ibam $bam -g $bed | grep 'chrM' - > $coverage
fi

#Maximum value of genome coverage in BED file
maximum=$(sort -nk 4 $coverage | tail -1 - | awk '{print $4}' -)

#Number of positions with X number of rNMPs
for i in $(seq 1 $maximum); do
	positions1+=($(awk '$4 == ('$i')' $coverage | wc -l))
done

#Number of positions with 0 rNMPs (total positions-positions with rNMPs)
positions2=$(echo "($total-$(wc -l $coverage | awk '{print $1}' -))" | bc)

#Print observed count data to fit to Poisson distribution
( IFS=$'\n'; echo -e "$positions2\n${positions1[*]}" ) > $counts1

lambda=$(echo "scale = 12; {positions1[*]}/2" | bc | awk '{printf "%.12f\n", $0}')
echo $lambda
