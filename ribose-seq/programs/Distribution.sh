#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program creates bedgraph files for forward and reverse strands

#Usage statement
function usage () {
	echo "Usage: Distribution.sh [options]
		-s Sample name(s) (e.g., FS1 FS2 FS3)
		-r Reference genome/Basename of Bowtie2 index
		-d Ribose-Map directory (e.g., /path/to/Ribose-Map)"
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
	#Input files
	bed=$directory/References/$reference.bed
	bam=$directory/Results/$reference/$sample/Alignment/$sample.bam
	coordinates=$directory/Results/$reference/$sample/Coordinates/$sample-Coordinates.bed
	
	#Output directory
	output=$directory/Results/$reference/$sample/Distribution
	
#############################################################################################################################
	#Create directory
	mkdir -p $output
	
	#Remove any old files
	rm -f $output/*.{bed,bg}

#############################################################################################################################
	if [[ -s $coordinates ]]; then
		
		#Count number of unique lines
		uniq -c $coordinates > $output/temp1.txt
		
		#Save coverage of rNMPs per chromosome to separate files
		for chr in $( awk '{print $1}' $bed ); do
			grep -w "$chr" $output/temp1.txt > $output/$sample-Distribution.$chr.bed
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
				
		#Remove temporary files
		rm -f $output/temp1.txt
		
		#Print completion status for program
		echo "Status: Program complete for $sample"
	
	fi
done
