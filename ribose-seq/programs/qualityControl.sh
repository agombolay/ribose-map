#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#############################################################################################################################
. /data2/users/agombolay3/Ribose-Map/config.txt
output=$directory/$name/pre-processing; mkdir -p $output

#Single-end reads
if [[ ! $read2 ]]; then
	fastqc $read1_fastq -o $directory/$name/pre-processing
	
	cutadapt $read1_fastq -m 50 -a $adapter -o $output/${name}_trimmed1.fq

#Paired-end reads
elif [[ $read2 ]]; then
	fastqc $read1_fastq $read2_fastq -o $directory/$name/pre-processing
	
	cutadapt $read1_fastq $read2_fastq -m 50 -a $adapter -o $output/${name}_trimmed1.fq \
	-p $output/${name}_trimmed1.fq
fi
