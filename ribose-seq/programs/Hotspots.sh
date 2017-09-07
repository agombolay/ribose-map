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
	#Input file
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-Coordinates.bed
	
	#Output directory
	output=$directory/Ribose-Map/Results/$reference/$sample/Hotspots
	
	#Create directory
	mkdir -p $output

#############################################################################################################################
	if [[ -s $coordinates ]]; then
		
		#Remove old files
		rm -f $output/$sample-*.bedgraph
		
		#Count number of unique lines
		uniq -c $coordinates > $output/temp1.txt

		#Add trackline for forward strand to input into UCSC genome browser
		echo "track type=bedGraph name="$sample-ForwardStrand" description="$sample-ForwardStrand" \
		color=0,128,0 visibility=full" > $output/$sample-Forward.bedgraph
		
		#Add trackline for reverse strand to input into UCSC genome browser
		echo "track type=bedGraph name="$sample-ReverseStrand" description="$sample-ReverseStrand" \
		color=0,128,0 visibility=full" > $output/$sample-Reverse.bedgraph
		
		#Rearrange file so format is same as bedgraph format (forward)
		awk -v "OFS=\t" '$5 == "+" {print $2, $3, $4, $1}' $output/temp1.txt >> $output/$sample-Forward.bedgraph

		#Rearrange file so format is same as bedgraph format (reverse)
		awk -v "OFS=\t" '$5 == "-" {print $2, $3, $4, $1}' $output/temp1.txt >> $output/$sample-Reverse.bedgraph

#############################################################################################################################
		#Remove old files
		rm -f $output/$sample-Coverage.bed $output/$sample-Hotspots.$chr.bed
		
		#Create file containing coverage of both strands
		awk -v "OFS=\t" '{print $2, $3, $4, $1}' temp1.txt > $output/$sample-Coverage.bed
		
		#Save coverage of rNMPs for each chromosome to separate files for plotting
		for chr in $( awk '{print $1}' $directory/Ribose-Map/References/$reference.bed ); do
			grep -w "$chr" $output/$sample-Coverage.bed > $output/$sample-Hotspots.$chr.bed
		done
		
#############################################################################################################################
		#Print completion status
		echo "Hotspots in $sample have been located"
		
		rm -f temp1.txt
	
	fi
done
