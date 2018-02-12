#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#############################################################################################################################
#Input files
fastq1=$directory/fastqs/$read1
fastq2=$directory/fastqs/$read2

#Remove old files
rm -f $output/*_trimmed{1,2}.fq

#############################################################################################################################
umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/extract1.fq

umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/extract1.fq --read2-in=$fastq2 --read2-out=$output/extract2.fq

grep -B 1 -A 2 ^$bc $output/extract1.fq | sed '/^--$/d' | awk 'NR%2 == 0 {sub(/^.{'${#bc}'}/,"")} {print}' > $output/filter.fq
