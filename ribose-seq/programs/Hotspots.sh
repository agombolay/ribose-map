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
	bed=$directory/Ribose-Map/References/$reference.bed
	bam=$directory/Ribose-Map/Results/$reference/$sample/Alignment/$sample.bam
	
	#Output files
	coverage=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Coverage.bed
	forward=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Forward.bedgraph
	reverse=$directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-Reverse.bedgraph
	
  #############################################################################################################################
	if [[ -s $bam ]]; then

		#Determine coverage at 3' position of reads
		bedtools genomecov -ibam $bam -d -3 > $coverage
  
		#Separate BAM file by forward and reverse strands
		samtools view -bS -f 16 $bam > reverse.bam
		samtools view -bS -F 16 $bam > forward.bam
		
		bedtools genomecov -bg -3 -trackline -trackopts 'name="ReverseStrand" description="Ribose-seq (Reverse)" \
		color=0,0,255 visibility=2' -ibam reverse.bam > $reverse
		
		bedtools genomecov -bg -3 -trackline -trackopts 'name="ForwardStrand" description="Ribose-seq (Forward)" \
		color=0,128,0 visibility=2' -ibam forward.bam > $forward
  
		#Save coverage of rNMPs per chromosome
		for chr in $( awk '{print $1}' $bed ); do
			grep -w "$chr" $coverage > $directory/Ribose-Map/Results/$reference/$sample/Hotspots/$sample-$chr.bed
		done

#############################################################################################################################
	#Print completion status
	echo "Hotspots in $sample have been located"
		
	rm -f reverse.bam forward.bam
	
	fi
done
