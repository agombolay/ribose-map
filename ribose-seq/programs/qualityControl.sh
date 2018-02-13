#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#############################################################################################################################
adapter='AGTTGCGACACGGATCTCTCA'

. /data2/users/agombolay3/Ribose-Map/config.txt

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
