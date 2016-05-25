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

output=$outputDirectory/transcribedRegions
if [[ ! -d $output ]];
then
	mkdir -p $output
fi

genes="$HOME/data/ribose-seq/results/transcribedRegions/genes.bed"
genes_Nuclear="$HOME/data/ribose-seq/results/transcribedRegions/genes.nuclear.bed"
genes_Mitochondrial="$HOME/data/ribose-seq/results/transcribedRegions/genes.mitochondrial.bed"

complementRegions="$HOME/data/ribose-seq/results/transcribedRegions/complementRegions.bed"
complementRegions_Nuclear="$HOME/data/ribose-seq/results/transcribedRegions/complementRegions.nuclear.bed"
complementRegions_Mitochondrial="$HOME/data/ribose-seq/results/transcribedRegions/complementRegions.mitochondrial.bed"

#TRANSCRIPTION ANALYSIS

#Sorts and merges overlapping genomic regions of sgdGene.bed file

for strand in $strands; do

	#Assign +/- value to variable, "strand"
	#Command accepts input from sgdGene.bed file
	awk -v strand=$strand '$6 == strand' < $sgdGene |
	
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

#Obtain genes in chromosome M
grep '^chrM' $genes > $genes_Mitochondrial

#Obtain non-coding genomic intervals in chromosome M
grep '^chrM' $complementRegions > $complementRegions_Mitochondrial

#Obtain genes in chromosomes I-XVI (not chromosome M or 2micron)
grep -v '^chrM' $genes | grep -v '^2micron' > "$genes_Nuclear"

#Obtain non-coding genomic intervals in chromosomes I-XVI (not	chromosome M or	2micron)
grep -v '^chrM' $complementRegions | grep -v '^2micron' > "$complementRegions_Nuclear"

echo "Determination of gene coordinates and complementary regions in genome complete"
