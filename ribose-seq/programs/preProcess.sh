#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

output=$directory/results/alignment

if [[ $fastq1 ]] && [[ ! $fastq2 ]]; then
  umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/umi_extracted1.fq
  
elif [[ $fastq1 ]] && [[ $fastq2 ]]; then
  umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/umi_extracted1.fq \
  --read2-in=$fastq2 --read2-out=$output/umi_extracted2.fq
fi

if [[ $barcode ]]; then
  #Filter reads by barcode and remove barcode from reads
  grep -B 1 -A 2 ^$barcode $output/umi_extracted1.fq | sed '/^--$/d' \
  | awk 'NR % 2 == 0 {sub(/^.{'${#barcode}'}/,"")} {print}' > $output/barcode.fq
  
  #Calculate % of reads that contain correct barcode sequence
  x=$(echo $(bc -l <<< "$(wc -l < $output/umi_extracted1.fq)/4")/$(bc -l <<< "$(wc -l < $output/barcode.fq)/4"))
  
  #Save info about % of reads that contain correct barcode sequence
  echo -e "Reads with the molecular barcode, $barcode: $(bc -l <<< "scale = 2; $x * 100")%" > $output/barcode.log
fi
