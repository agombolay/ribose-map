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
	
	for subset in "mito" "nucleus"; do
		
		#Directory
		mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Distribution
	
		#Input file
		bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample-MappedReads.bam
	
		#Output file
		counts=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Counts.$subset.txt
		coverage=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Coverage.$subset.bed
	
		#Remove old files
		rm -f $coverage $counts temp{1..3}.txt
	
		#Determine coverage at 3' position of reads
		bedtools genomecov -d -3 -ibam $bam > temp1.txt
	
		#Select region of genome (i.e., nucleus or mito)
		if [ $subset == "mito" ]; then
			grep -E '(chrM|MT)' temp1.txt > $coverage
		elif [ $subset == "nucleus" ]; then
			grep -vE '(chrM|MT)' temp1.txt > $coverage
		fi
		
		#Sort by # of rNMPs
		sort -k3n $coverage > temp2.txt

		#Maximum # of rNMPs in observed data
		max=$(tail -1 temp2.txt | awk '{print $3}' -)

		#Number of positions containing 0...max # of rNMPs
		for i in $(seq 0 $max); do
			awk '$3 == ('$i')' temp2.txt | wc -l >> temp3.txt
		done
		
#############################################################################################################################
		#STEP 2: Create dataset file of observed rNMP counts
	
		#Add column names to header line
		echo -e "rNMPs\tPositions" > $counts

		#Add number of positions containing 0...max # of rNMPs
		paste <(echo "$(seq 0 $max)") <(cat temp3.txt) >> $counts

		#Print completion status
		echo "Counts for $sample ($subset) have been determined"
	
		#Remove temp files
		rm -f temp{1..3}.txt

	done
done
