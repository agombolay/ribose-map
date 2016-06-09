#!/usr/bin/env bash

#Author: Alli Gombolay 
#Script to download the reference genome files required for the Ribose-seq Analysis Pipeline

#Download UCSC's twoBitToFa program to convert .2bit file to .fa
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

#Download .2bit file of the reference genomic sequence from UCSC's site
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.2bit

#Convert the reference genome sequence file from .2bit to .fa
twoBitToFa sacCer2.2bit sacCer2.fa

#Build Bowtie index for the reference genome from the .fa file
bowtie-build sacCer2.fa sacCer2Index

#Download file of size (in base pairs) of the reference genome from UCSC's site
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.chrom.sizes

#Sort the reference genome file for processing
sort sacCer2.chrom.sizes -o sacCer2.chrom.sizes

#OR

#Download fetchChromSizes to create file containing sizes of chromosomes
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
#fetchChromSizes sacCer2 > sacCer2.chrom.sizes

#Download file containing locations of genes on chromosomes from UCSC's site
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/database/sgdGene.txt.gz

#Uncompress the .txt.gz file and then convert it to .bed file format (rearrange columns and remove some)
gunzip sgdGene.txt.gz | cat sgdGene.txt | awk  ' {OFS="\t"; print $3,$5,$6,".", ".",$4 } ' > sgdGene.bed
