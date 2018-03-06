#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#############################################################################################################################
#Load config file
. "$1"

#Create output directory
output=$directory/trimmed; mkdir -p $output

#############################################################################################################################

if [[ ! $instrument ]]; then
	nextseq=''
elif [[ $instrument ]]; then
	nextseq='--nextseq-trim=20'
fi

#Single-end reads
if [[ ! $raw2 ]]; then
	fastqc "$2" -o $output
	cutadapt $nextseq -a $adapter -m 50 "$2" -q 10 -o $output/$sample-trimmed.fastq
	fastqc $output/$sample-trimmed.fastq -o $output

#Paired-end reads
elif [[ $read2 ]]; then
	fastqc "$2" "$3" -o $output
	cutadapt $nextseq -a $adapter -m 50 "$2" "$3" -q 10 -o $output/$sample-trimmed1.fastq -p $output/$sample-trimmed2.fastq
	fastqc $output/$sample-trimmed1.fastq $output/$sample-trimmed2.fastq -o $output

fi
