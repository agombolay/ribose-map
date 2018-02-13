#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

. /data2/users/agombolay3/Ribose-Map/config.txt

#############################################################################################################################
#Input files
#fastq1=$directory/fastqs/$read1
#fastq2=$directory/fastqs/$read2

#############################################################################################################################
if [[ ! $read2 ]]; then
	#Single-end reads
	fastqc $fastq1 -o $directory/$name/alignment
	#cutadapt $fastq1 -m 50 -a 'AGTTGCGACACGGATCTCTCA' -o $directory/fastqs/${sample}_trimmed1.fq

elif [[ $read2 ]]; then
	#Paired-end reads
	fastqc $fastq1 $fastq2 -o $output
	cutadapt $fastq1 $fastq2 -m 50 -a 'AGTTGCGACACGGATCTCTCA' -o $directory/fastqs/$sample_trimmed1.fq -p $directory/fastqs/$sample_trimmed2.fq
fi
