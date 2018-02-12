#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#############################################################################################################################
$output=$directory/results/alignment

if [[ $fastq1 ]]; then
umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/UMI1.fq
umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/UMI1.fq --read2-in=$fastq2 --read2-out=$output/UMI2.fq
grep -B 1 -A 2 ^$bc $output/UMI1.fq | sed '/^--$/d' | awk 'NR%2==0 {sub(/^.{'${#bc}'}/,"")} {print}' > $output/bc.fq
