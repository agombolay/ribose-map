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
	
	for subset in "2micron" "nucleus"; do
		
		#Create directory
		mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Distribution
	
		#Input file
		bed=$directory/Ribose-Map/References/$reference.bed
		bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample.bam
		coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-Coordinates.$subset.bed
		
		#Output file
		windows=$directory/Ribose-Map/References/$reference-windows.bed
		coverage=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Coverage.bed
		
		#Select region of genome (i.e., nucleus or mito)
		if [ $subset == "2micron" ]; then
		
			bedtools makewindows -g $bed -w 100 > $windows
			bedtools intersect -a $windows -b $coordinates -c -nonamecheck > temp1.bed
		
			grep -E '(chrM|MT)' temp1.bed > temp2.bed
		
				#Sort by # of rNMPs
				sort -k4n temp2.bed > temp3.bed

				#Maximum # of rNMPs in observed data
				max=$(tail -1 temp3.bed | awk '{print $4}' -)

				#Number of positions containing 0...max # of rNMPs
				for i in $(seq 0 $max); do
					awk '$4 == ('$i')' temp3.bed | wc -l >> temp3.txt
				done
		
				#STEP 2: Create dataset file of observed rNMP counts
				counts=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Counts.$subset.txt
				#Add column names to header line
				echo -e "rNMPs\tWindows" > $counts

				#Add number of positions containing 0...max # of rNMPs
				paste <(echo "$(seq 0 $max)") <(cat temp3.txt) >> $counts

				#Print completion status
				echo "Counts for $sample ($chr) have been determined"
	
				#Remove temp files
				rm -f temp3.txt
				
		elif [ $subset == "nucleus" ]; then
			
			bedtools makewindows -g $bed -w 2500 > $windows
			bedtools intersect -a $windows -b $coordinates -c -nonamecheck > temp1.bed
		
			for chr in $( awk '{print $1}' $bed ); do
				
				counts=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Counts.$chr.txt
				grep -w "$chr" temp1.bed > temp2.bed
		
				#Sort by # of rNMPs
				sort -k4n temp2.bed > temp3.bed

				#Maximum # of rNMPs in observed data
				max=$(tail -1 temp3.bed | awk '{print $4}' -)

				#Number of positions containing 0...max # of rNMPs
				for i in $(seq 0 $max); do
					awk '$4 == ('$i')' temp3.bed | wc -l >> temp3.txt
				done
		
#############################################################################################################################
				#STEP 2: Create dataset file of observed rNMP counts
	
				#Add column names to header line
				echo -e "rNMPs\tWindows" > $counts

				#Add number of positions containing 0...max # of rNMPs
				paste <(echo "$(seq 0 $max)") <(cat temp3.txt) >> $counts

				#Print completion status
				echo "Counts for $sample ($chr) have been determined"
	
				#Remove temp files
				rm -f temp3.txt
		
			done
		fi
	done
done
