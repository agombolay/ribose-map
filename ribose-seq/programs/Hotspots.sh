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

	#Create directory
	mkdir -p $directory/Ribose-Map/Results/$reference/$sample/Hotspots

	#Input files
	bed=$directory/Ribose-Map/References/$reference.bed
	bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample.bam
	coverage=$directory/Ribose-Map/Results/$reference/$sample/Distribution/$sample-Coverage.bed
	
	if [[ -s $coverage ]] && [[ -s $bam ]]; then
		
		#Output files
		forward=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Forward.bedgraph
		reverse=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Reverse.bedgraph
	
		#Determine coverage at 3' position of reads
		samtools view -f 16 $bam | bedtools genomecov -bg -3 -trackline -trackopts color=0,0,255 -ibam - > $forward
		samtools view -F 16 $bam | bedtools genomecov -bg -3 -trackline -trackopts color=0,128,0 -ibam - > $reverse

		#Select chromosomes
		genome=$(awk '{print $1}' $bed)	
	
		#Save coverage of rNMPs to separate files per chromosome
        	for chr in ${genome[@]}; do
			hotspots=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Hotspots.$chr.bed
			grep -w "$chr" $coverage > $directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-$chr.bed
		done
	
		#Print completion status
		echo "Hotspots in $sample have been located"
	
	fi
done
