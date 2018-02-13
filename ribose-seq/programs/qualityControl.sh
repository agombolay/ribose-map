#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#############################################################################################################################
. /data2/users/agombolay3/Ribose-Map/config.txt

output=$directory/results/$sample/pre-processing
mkdir -p $output

#Single-end reads
if [[ ! $read2 ]]; then
	trim_galore --fastqc -a 'AGTTGCGACACGGATCTCTCA' --length 50 $read1_fastq -o $output

#Paired-end reads
elif [[ $read2 ]]; then
	trim_galore --fastqc -a 'AGTTGCGACACGGATCTCTCA' --length 50 --paired $read1_fastq $read2_fastq -o $output

fi
