#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: July 11, 2016
#This program determines 0-based and 1-based coordinates of rNMPs (5’ end of each mapped read) for positive and negative strands

#Coordinates (0-based) of Sequencing Reads

#Covert file from BAM to BED format using BEDtools software
bedtools bamtobed -i data.bam > data.bed

#Convert file from BAM to SAM format using SAMtools software
samtools view data.bam > data.sam

#Extract coordinate information and save it to a new file, “data.coordinates.bed”
paste data.bed data.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' > data.coordinates.bed

#0-BASED COORDINATES OF rNMPs

#Obtain positions of rNMPs (5’ end of each mapped read) for positive strand:
bedtools genomecov -5 -strand + -bg -ibam FS15_processed.bam > FS15.rNMPs.positive.txt

#Obtain positions of rNMPs (5’ end of each mapped read) for negative strand:
bedtools genomecov -5 -strand - -bg -ibam FS15_processed.bam > FS15.rNMPs.negative.txt

#1-BASED COORDINATES OF	rNMPs

#Obtain positions of rNMPs (5’ end of each mapped read) for positive strand:
bedtools genomecov -5 -strand + -d -ibam FS15_processed.bam > FS15.rNMPs.positive.txt

#Remove rows where genome coverage equals 0
awk '$3 != 0' FS15.rNMPs.positive.txt > temporary

#Change filename back to original
mv temporary FS15.rNMPs.positive.txt

#Obtain positions of rNMPs (5’ end of each mapped read) for negative strand:
bedtools genomecov -5 -strand - -d -ibam FS15_processed.bam > FS15.rNMPs.negative.txt

#Remove rows where genome coverage equals 0
awk '$3 != 0' FS15.rNMPs.negative.txt > temporary

#Change filename back to original
mv temporary FS15.rNMPs.negative.txt



