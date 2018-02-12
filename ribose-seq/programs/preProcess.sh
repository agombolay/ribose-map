#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#############################################################################################################################
$output=$directory/results/alignment

if [[ $fastq1 ]] && [[ ! $fastq2 ]]; then
  umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/UMI1.fq
  
  if [[ $barcode ]]; then
  
    #Filter reads by barcode and remove barcode from reads
    grep -B 1 -A 2 ^$bc $output/UMI1.fq | sed '/^--$/d' | awk 'NR%2==0 {sub(/^.{'${#bc}'}/,"")} {print}' > $output/bc.fq
    
    #Calculate % of reads that contain correct barcode sequence
    x=$(echo $((`wc -l < $output/filtered.fq` / 4))/$((`wc -l < $output/UMI.fq` / 4)))

    #Save info about % of reads that contain correct barcode sequence
    echo -e "Percentage of reads with barcode: $(echo "$x*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/barcode.log

  fi
  
elif [[ $fastq1 ]] && [[ $fastq2 ]]; then
  umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/UMI1.fq --read2-in=$fastq2 --read2-out=$output/UMI2.fq
  
  if [[ $barcode ]]; then
  
    #Filter reads by barcode and remove barcode from reads
    grep -B 1 -A 2 ^$bc $output/UMI1.fq | sed '/^--$/d' | awk 'NR%2==0 {sub(/^.{'${#bc}'}/,"")} {print}' > $output/bc.fq
    
    #Calculate % of reads that contain correct barcode sequence
    x=$(echo $((`wc -l < $output/filtered.fq` / 4))/$((`wc -l < $output/UMI.fq` / 4)))

    #Save info about % of reads that contain correct barcode sequence
    echo -e "Percentage of reads with barcode: $(echo "$x*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/barcode.log

  fi
 
fi
