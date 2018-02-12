#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#Single-end reads
trim_galore $fastq1 --fastqc --length $minimum -a $adapter -o $output

#Paired-end reads
trim_galore $fastq1 $fastq2 --fastqc --paired --length $minimum -a $adapter $output/filtered.fq -o $output
