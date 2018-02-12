#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#Usage statement
function usage () {
	echo "Usage: pre-process.sh [options]
		Required:
		-d Ribose-Map repository
		-s Name of sequenced library
		-f Input Read 1 FASTQ filename
		-r Input Read 2 FASTQ filename"
}

while getopts "f:r:s:d:h" opt; do
    	case "$opt" in
		f ) read1=$OPTARG ;;
		r ) read2=$OPTARG ;;
		s ) sample=$OPTARG ;;
		d ) directory=$OPTARG ;;
    #Print usage statement
    h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Single-end reads
trim_galore $fastq1 --fastqc --length $minimum -a $adapter -o $output

#Paired-end reads
trim_galore $fastq1 $fastq2 --fastqc --paired --length $minimum -a $adapter $output/filtered.fq -o $output
