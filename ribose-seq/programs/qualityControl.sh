#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#Usage statement
function usage () {
	echo "Usage: qualityControl.sh [options]
		-d Ribose-Map directory
		-s Name of sequenced library
		-n Sequencing instrument used"
}

#Command-line options
while getopts "s:d:i:h" opt; do
    case "$opt" in
        s ) sample=$OPTARG ;;
	d ) directory=$OPTARG ;;
	i ) instrument=$OPTARG ;;
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
if [[ ! $nextseq ]]: then

	if [[ ! $read2 ]]; then
		fastqc $read1_fastq -o $output
		cutadapt -a $adapter -m 50 $read1_fastq -o $output/trimmed.fq

	#Paired-end reads
	elif [[ $read2 ]]; then
		fastqc $read1 $read2 -o $output
		cutadapt -a $adapter -m 50 -p $read1 $read2 -o $output/trimmed1.fq -p $output/trimmed2.fq

	fi
	
elif [[ ! $nextseq ]]: then

	if [[ ! $read2 ]]; then
		fastqc $read1_fastq -o $output
		cutadapt --nextseq-trim=20 -a $adapter -m 50 $read1_fastq -o $output/trimmed.fq

	#Paired-end reads
	elif [[ $read2 ]]; then
		fastqc $read1 $read2 -o $output
		cutadapt --nextseq-trim=20 -a $adapter -m 50 -p $read1 $read2 -o $output/trimmed1.fq -p $output/trimmed2.fq
	fi
fi
