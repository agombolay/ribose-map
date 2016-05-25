#! /usr/bin/env bash

#Author: Alli Gombolay
#This script determines the gene coordinates and non-coding genomic intervals in sacCer2 genome
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program (1_alignment.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-r] 'reference genome' [-d] 'Ribose-seq directory' [-h]
          -r Reference genome of interest (i.e., sacCer2)
          -d Location to save local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-i], [-d], and [-h])
while getopts "i:d:h" opt;
do
    case $opt in
	#Specify input as variable to allow only one input argument
        r ) reference=$OPTARG ;;
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

#VARIABLE SPECIFICATION

#DNA strands, positive or negative
strands="+ -"

#INPUT
#Input files for this script: sgdGene.bed and sacCer2.chrom.sizes from genome.ucsc.edu/index.html
#Please see $HOME/data/ribose-seq/data/reference/downloadReferenceGenome.sh for more information

#Location of file containing gene coordinates 
BED=$directory/reference/$reference.bed

#Location of file containing chromosome sizes
sizes=$directory/reference/$reference.chrom.sizes

#OUTPUT
#Location of output "ribose-seq" alignment directory
output=$directory/ribose-seq/results/transcribedRegions

#Create directory for output if it does not already exist
if [[ ! -d $output ]];
then
	mkdir -p $output
fi

#Location of output files containing genes
genes="$output/$(basename $BED .bed).nuclear.bed"

#Location of output files containing complement regions
complement="$output/$(basename $BED .bed).complement.bed"

#TRANSCRIPTION ANALYSIS

#Sorts and merges overlapping genomic regions of sgdGene.bed file
for strand in $strands; do

	#Assign +/- value to variable, "strand"
	#Command accepts input from sgdGene.bed file
	awk -v strand=$strand '$6 == strand' < $BED |
	
		#Takes output from previous step as input and sorts it alphanumerically
		bedtools sort -i - |

		#Takes output from previous step as input and merges overlapping genes
		bedtools merge -i - |
		
		#Prints chromosme names, start and end positions of genes, and strand
		awk -v strand=$strand 'BEGIN {OFS="\t"} {print $0, ".", ".", strand}'

#Takes output from previous step as input and executes a final sort
done | bedtools sort -i - > $genes

#Returns all genomic intervals not covered by at least one interval in input file
#Next, sorts output and prints chromosome names, start and end positions of genes
bedtools complement -i $genes -g $sizes | bedtools sort -i - |
awk 'BEGIN {OFS="\t"} {print $0, ".", ".", "."}' > $complement

#Obtain genes on chromosome M
grep '^chrM' $genes > "$output/$(basename $genes .bed).mito.bed"

#Obtain non-coding genomic intervals in chromosome M
grep '^chrM' $complement > "$output/$(basename $complement .bed).mito.bed"

#Obtain genes on chromosomes I-XVI (not on chromosome M or 2micron)
grep -v '^chrM' $genes | grep -v '^2micron' > "$output/$(basename $genes .bed).nuc.bed"

#Obtain non-coding genomic intervals in chromosomes I-XVI (not on chromosome M or 2micron)
grep -v '^chrM' $complement | grep -v '^2micron' > "$output/$(basename $complement .bed).nuc.bed"
