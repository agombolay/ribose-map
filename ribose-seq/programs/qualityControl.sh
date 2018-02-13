#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#Usage statement
function usage () {
	echo "Usage: qualityControl.sh [options]
		-f Filepath of read 1
		-r Filepath of read 2
		-d Ribose-Map directory
		-s Name of sequenced library
		-n Sequencing instrument used"
}

#Command-line options
while getopts "s:f:r:d:i:h" opt; do
    case "$opt" in
	s ) sample=$OPTARG ;;
	f ) forward=$OPTARG ;;
	r ) reverse=$OPTARG ;;
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
output=$directory/results/$sample/alignment
mkdir -p $output

adapter='AGTTGCGACACGGATCTCTCA'
#############################################################################################################################

if [[ ! $instrument ]]; then
	nextseq=''
elif [[ $instrument ]]; then
	nextseq='--nextseq-trim=20'
fi

#Single-end reads
if [[ ! $read2 ]]; then
	fastqc $forward -o $output
	
	cutadapt $nextseq -a $adapter -m 50 $forward -o $output/qc.fq
	
	fastqc $output/qc.fq -o $output

#Paired-end reads
elif [[ $read2 ]]; then
	fastqc $forward $reverse -o $output
	
	cutadapt $nextseq -a $adapter -m 50 $forward $reverse -o $output/qc1.fq -p $output/qc2.fq
	
	fastqc $output/qc1.fq -$output/qc2.fq -o $output

fi
