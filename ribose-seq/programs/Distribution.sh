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

#############################################################################################################################
	
	#STEP 1: Count number of positions containing 0...max # of rNMPs
	
	#Directory
	mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Distribution
	
	#Input file
	bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample-MappedReads.bam
		
	#Output file
	coverage=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Coverage.bed

	#Remove old files
	rm -f $coverage temp{1..2}.txt
		
	#Determine coverage at 3' position of reads
	bedtools genomecov -d -3 -ibam $bam > temp1.txt
	
	for subset in "mito" "nucleus"; do
		
		#Output file
		counts=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Counts.$subset.txt
	
		#Remove old file
		rm -f $counts
	
		#Select regions of interest and sort by # of rNMPs
		if [ $subset == "mito" ]; then
			grep -E '(chrM|MT)' $coverage | sort -k3n - > temp1.txt
		elif [ $subset == "nucleus" ]; then
			grep -vE '(chrM|MT)' $coverage | sort -k3n - > temp1.txt
		fi

		#Maximum # of rNMPs in observed data
		max=$(tail -1 temp1.txt | awk '{print $3}' -)

		#Number of positions containing 0...max # of rNMPs
		for i in $(seq 0 $max); do
			awk '$3 == ('$i')' temp1.txt | wc -l >> temp2.txt
		done
		
#############################################################################################################################
		#STEP 2: Create dataset file of observed rNMP counts
	
		#Add column names to header line
		echo -e "rNMPs\tPositions" > $counts

		#Add number of positions containing 0...max # of rNMPs
		paste <(echo "$(seq 0 $max)") <(cat temp2.txt) >> $counts

		#Print completion status
		echo "Counts for $sample ($subset) have been determined"
	
		#Remove temp files
		rm -f temp{1..2}.txt

	done
done
