#!/usr/bin/env bash

#Author: Alli Gombolay 
#Script to download the reference genome files required for the Ribose-seq Analysis Pipeline

#Download UCSC's twoBitToFa program to convert .2bit files to .fa
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

#Download the .2bit file of the reference genome sequence from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.2bit

#Convert the reference genome sequence file from .2bit to .fa
twoBitToFa sacCer2.2bit sacCer2.fa

#Build Bowtie index for the reference genome from the .fa file
bowtie-build sacCer2.fa sacCer2Index

#Download file containing length in base pairs of the reference genome from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.chrom.sizes

#Sort the reference genome file for processing
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
