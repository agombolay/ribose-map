#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

path="$1"
. $path/config.txt

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
	fastqc "$2" -o $output
	cutadapt $nextseq -a $adapter -m 50 "$2" -o $output/$sample-trimmed.fastq
	fastqc $output/$sample-trimmed.fq -o $output

#Paired-end reads
elif [[ $read2 ]]; then
	fastqc "$2" "$3" -o $output
	cutadapt $nextseq -a $adapter -m 50 "$2" "$3" -o $output/$sample-trimmed1.fastq -p $output/$sample-trimmed2.fastq
	fastqc $output/$sample-trimmed1.fq $output/$sample-trimmed2.fq -o $output

fi
