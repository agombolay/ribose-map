#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

. /data2/users/agombolay3/Ribose-Map/config.txt

output=$directory/results/$sample/alignment
mkdir -p $output

#############################################################################################################################

if [[ ! $instrument ]]; then
	nextseq=''
elif [[ $instrument ]]; then
	nextseq='--nextseq-trim=20'
fi

#Single-end reads
if [[ ! $raw2 ]]; then
	fastqc $read1 -o $output
	cutadapt $nextseq -a $adapter -m 50 $raw1 -o $output/${sample}_trimmed.fq
	fastqc $output/${sample}_trimmed.fq -o $output

#Paired-end reads
elif [[ $read2 ]]; then
	fastqc $read1 $read2 -o $output
	cutadapt $nextseq -a $adapter -m 50 $raw1 $raw2 -o $output/${sample}_trimmed1.fq -p $output/${sample}_trimmed2.fq
	fastqc $output/${sample}_trimmed1.fq $output/${sample}_trimmed2.fq -o $output

fi
