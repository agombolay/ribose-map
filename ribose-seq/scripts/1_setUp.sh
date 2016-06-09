#!/usr/bin/env bash

#Author: Alli Gombolay 
#Script to download the reference genome files required for the Ribose-seq Analysis Pipeline

#COMMAND LINE OPTIONS

#Name of the program (1_setUp.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-r] 'reference genome' [-d] 'directory' [-h]
          -r referene genome of interest
          -d"

}

#Use getopts function to create the command-line options ([-i], [-b], [-d], and [-h])
while getopts "r:h" opt;
do
    case $opt in
	#Specify input as variable to allow only one input argument
	b ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

#Download twoBit file of the reference genome sequence from UCSC's site
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.2bit

#Convert the reference genome sequence file from .2bit to .fa
twoBitToFa sacCer2.2bit sacCer2.fa

#Build Bowtie index for the reference genome from the .fa file
bowtie-build sacCer2.fa sacCer2Index

#Download file of chromosome sizes (bp) of the reference genome from UCSC's site
#Note: fetchChromSizes is a UCSC program that can also be used to create the file
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.chrom.sizes

#Sort the reference genome file for processing
sort sacCer2.chrom.sizes -o sacCer2.chrom.sizes

#Download file of gene locations (start and end positions) from UCSC's site
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/database/sgdGene.txt.gz

#Uncompress the .txt.gz file and then convert it from .txt to .bed (rearrange columns and remove some)
gunzip sgdGene.txt.gz | cat sgdGene.txt | awk  ' {OFS="\t"; print $3,$5,$6,".", ".",$4 } ' > sgdGene.bed
