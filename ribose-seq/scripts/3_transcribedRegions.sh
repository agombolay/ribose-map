#! /usr/bin/env bash

#Author: Alli Gombolay
#This script determines the gene coordinates and non-coding genomic intervals in sacCer2 genome
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#VARIABLE SPECIFICATION

#DNA strands, positive or negative
strands="+ -"

#INPUT
#Input files for this script: sgdGene.bed and sacCer2.chrom.sizes from genome.ucsc.edu/index.html
#Please see $HOME/data/ribose-seq/data/reference/downloadReferenceGenome.sh for more information

#Location of sgdGene.bed file
sgdGene=$HOME/data/ribose-seq/data/reference/sgdGene.bed

#Location of sacCer2.chrom.sizes file
chromosomeSizes=$HOME/data/ribose-seq/data/reference/sacCer2.chrom.sizes

#OUTPUT
#Location of output "ribose-seq" alignment directory
output=$directory/ribose-seq/results/transcribedRegions

#Create directory for output
if [[ ! -d $output ]];
then
	mkdir -p $output
fi

BED=$directory/genes.bed
genes="$output/$(basename $genesbed .bed).nuclear.bed"
$complementRegions="$output/$(basename $genesbed .bed).complementRegions.bed"

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
bedtools complement -i $genes -g $chromosomeSizes | bedtools sort -i - |
awk 'BEGIN {OFS="\t"} {print $0, ".", ".", "."}' > $complementRegions

#Obtain genes on chromosome M
grep '^chrM' $genes > "$output/$(basename $genes .bed).mito.bed"

#Obtain non-coding genomic intervals in chromosome M
grep '^chrM' $complementRegions > "$output/$(basename $complementRegions .bed).mito.bed"

#Obtain genes on chromosomes I-XVI (not chromosome M or 2micron)
grep -v '^chrM' $genes | grep -v '^2micron' > "$output/$(basename $genes .bed).nuc.bed"

#Obtain non-coding genomic intervals in chromosomes I-XVI (not	chromosome M or	2micron)
grep -v '^chrM' $complementRegions | grep -v '^2micron' > "$output/$(basename $complementRegions .bed).nuc.bed"

echo "Determination of gene coordinates and complementary regions in genome complete"
