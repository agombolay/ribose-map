#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates the coveage at each rNMP position

#Usage statement
function usage () {
	echo "Usage: Hotspots.sh [options]
		-s Sample name(s) (e.g., FS1, FS2, FS3)
		-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:r:d:h" opt; do
    case "$opt" in
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
	#Create directory
	mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Hotspots
	
	#Input files
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-Coordinates.$subset.bed
	
	#Output files
	coverage=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Coverage.bed
	forward=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Forward.bedgraph
	reverse=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Reverse.bedgraph
	
  ###########################################################################################################################
	if [[ -s $coordinates ]]; then
		
		#Add trackline for forward strand to input into UCSC genome browser
		echo "track type=bedGraph name="ForwardStrand" description="$sample" color=0,128,0 visibility=2" > $forward
		
		#Create bedgraph file for reverse strand to input into UCSC genome browser
		echo "track type=bedGraph name="ForwardStrand" description="$sample" color=0,128,0 visibility=2" > $reverse
		
		#Count number of unique lines
		uniq -c $coordinates > temp1.txt
		
		#Rearrange file so format is same as bedgraph format (forward)
		awk -v "OFS=\t" '$5 == "+" {print $2, $3, $4, $1}' temp1.txt) >> $forward

		#Rearrange file so format is same as bedgraph format (reverse)
		awk -v "OFS=\t" '$5 == "-" {print $2, $3, $4, $1}' temp1.txt) >> $reverse
		
		#Save coverage of rNMPs per chromosome
		for chr in $( awk '{print $1}' $bed ); do
			hotspots=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Hotspots.$chr.bed
			grep -w "$chr" $coverage > $hotspots
		done
		
#############################################################################################################################
	#Print completion status
	echo "Hotspots in $sample have been located"
		
	rm -f temp1.txt
	
	fi
done
