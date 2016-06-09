#!/usr/bin/env bash

#Author: Alli Gombolay 
#Script to download reference genome (sacCer2) files

#Download twoBitToFa program to convert .2bit file to .fa
#This program is necessary to convert sacCer2.2bit to .fa format
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

#Download the sacCer2 reference genome (Chromosomes I-XVI, ChrM, and 2micron plasmid)
#sacCer2.2bit contains the sequences of all 16 chromosomes, chr M, and 2micron plasmid
wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/bigZips/sacCer2.2bit

#Run twoBitToFa program to convert file from twoBit to FASTA format
#sacCer2.2bit file contains the genomic sequence of the reference, sacCer2
twoBitToFa sacCer2.2bit sacCer2.fa

#Build Bowtie index for sacCer2
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
