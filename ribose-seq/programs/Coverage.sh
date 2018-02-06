#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Creates bedgraph files for forward and reverse strands
#2. Saves coverage of rNMPs per chromosome to separate files

#Usage statement
function usage () {
	echo "Usage: Coverage.sh [options]
		-d Ribose-Map repository
		-s Name of sequenced library
		-i Basename of Bowtie2 index"
}

#Command-line options
while getopts "d:s:i:h" opt; do
    case "$opt" in
        s ) sample=$OPTARG ;;
	i ) reference=$OPTARG ;;
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
reference=$directory/References/$reference.bed
coordinates=$directory/Results/$reference/$sample/Coordinates/$sample-Coordinates.bed
	
#Create output directory and remove old directory if present
output=$directory/Results/$reference/$sample/Coverage; mkdir -p $output; rm -rf $output

#############################################################################################################################
if [[ -s $coordinates ]]; then
		
	#Save coverage of rNMPs per chromosome to separate files
	for chr in $( awk '{print $1}' $bed ); do
		uniq -c $coordinates | grep -w "$chr" - > $output/$sample-Coverage.$chr.bed
	done
		
	#Add trackline for forward strand to input into UCSC genome browser
	echo "track type=bedGraph name="$sample-ForwardStrand" description="$sample-ForwardStrand" \
	color=0,128,0 visibility=full" > $output/$sample-Forward.bg
		
	#Add trackline for reverse strand to input into UCSC genome browser
	echo "track type=bedGraph name="$sample-ReverseStrand" description="$sample-ReverseStrand" \
	color=0,0,255 visibility=full" > $output/$sample-Reverse.bg
		
	#Rearrange file so format is same as bedgraph format (forward)
	awk -v "OFS=\t" '$5 == "+" {print $2,$3,$4,$1}' $output/temp1.txt >> $output/$sample-Forward.bg

	#Rearrange file so format is same as bedgraph format (reverse)
	awk -v "OFS=\t" '$5 == "-" {print $2,$3,$4,$1}' $output/temp1.txt >> $output/$sample-Reverse.bg

#############################################################################################################################
	#Print completion status
	echo "Status: Program complete for $sample"
	
fi
