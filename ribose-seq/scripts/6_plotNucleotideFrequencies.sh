#! /usr/bin/env bash

#Author: Alli Gombolay
#This program plots the frequencies of the nucleotides located in the assembled samples.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

offset_values=(100 50 15)

modes=("all" "only-mitochondria" "no-mitochondria" "only-2micron")

input=$directory/ribose-seq/results/$samples/nucleotideFrequencies/

output="$directory/ribose-seq/results/$samples/plots/nucleotideFrequencies"

if [[ ! -d $output ]];
then
  mkdir -p $output
fi

for index in ${!modes[@]};
do

  mode=${modes[$index]}

  for value in ${offset_values[@]};
  do
    
    sampleID="$sample.subset-$ignore_mode"
    
    tables="input/$sample.$mode.nucleotideFrequencies.tab"
    
    Rscript 6_nucleotideFrequencies.R -n "$sampleID" -d $output --offsetmax $value $tables
    
  done

done
