#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#Usage statement
function usage () {
	echo "Usage: qualityControl.sh [options]
		-d Ribose-Map directory
		-s Name of sequenced library
		-a Sequencing adapter to trim"
}

#Command-line options
while getopts "s:d:a:h" opt; do
    case "$opt" in
        s ) sample=$OPTARG ;;
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
	h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#############################################################################################################################
adapter='AGTTGCGACACGGATCTCTCA'

output=$directory/results/$sample/alignment
mkdir -p $output

#Single-end reads
if [[ ! $read2 ]]; then
	fastqc $read1_fastq -o $output
	cutadapt-a $adapter -m 50 $read1_fastq -o $output/${sample}_trimmed.fq

#Paired-end reads
elif [[ $read2 ]]; then
	fastqc $read1 $read2 -o $output
	cutadapt -a $adapter -m --paired $read1 $read2 -o $output/${sample}_trimmed1.fq -p $output/${sample}_trimmed2.fq

fi

--nextseq-trim=20
