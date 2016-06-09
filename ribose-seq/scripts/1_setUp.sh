#!/usr/bin/env bash

#Author: Alli Gombolay 
#Script to download reference genome (sacCer2) files

#Download twoBitToFa program to convert twoBit (.2bit) files to FASTA (.fa)
#This program is necessary to convert the reference genome sequence file to .fa
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

#Download the reference genome sequence .2bit file from the UCSC Genome Browser
#sacCer2.2bit contains the sequences of Chromosomes I-XVI, ChrM, and 2micron plasmid
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.2bit

#Run twoBitToFa program to convert reference genome sequence file from .2bit to .fa
twoBitToFa sacCer2.2bit sacCer2.fa

#Build the Bowtie index for the reference genome from the output .fa file
bowtie-build sacCer2.fa sacCer2Index

#Download file containing sizes of chromosomes of sacCer2 genome
#Contains length in base pairs of all 16 chromosomes, chr M, and 2micron plasmid
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.chrom.sizes

sort sacCer2.chrom.sizes -o sacCer2.chrom.sizes

#OR

#Download fetchChromSizes to create file containing sizes of chromosomes
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
#fetchChromSizes sacCer2 > sacCer2.chrom.sizes

#Download sgdGene.txt.gz and convert it to BED format
#Contains names of chromosomes, start and end positions of genes, and strand
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/database/sgdGene.txt.gz

#Uncompress sgdGene.txt.gz file and then convert it to a bed file (rearrange columns and remove some)
gunzip sgdGene.txt.gz | cat sgdGene.txt | awk  ' {OFS="\t"; print $3,$5,$6,".", ".",$4 } ' > sgdGene.bed
