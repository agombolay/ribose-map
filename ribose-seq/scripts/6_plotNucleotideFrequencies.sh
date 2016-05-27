#! /usr/bin/env bash

#Author: Alli Gombolay
#This program plots the frequencies of the nucleotides located in the assembled samples.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

offset_values=(100 50 15)

ignore_modes=("all" "only-mito" "no-mito" "only-2micron")

input=$directory/ribose-seq/results/$samples/nucleotideFrequencies/
output="$directory/ribose-seq/results/$samples/plots/nucleotideFrequencies"

if [[ ! -d $output ]];
then
  mkdir -p $output
fi

for index in ${!ignore_modes[@]};
do

  ignore_mode=${ignore_modes[$index]}

  for value in ${offset_values[@]};
  do
    sampleID="$sample.subset-$ignore_mode"
    tables="input/$sample.ignore.$ignore_mode.nucleotideFrequencies.tab"
    Rscript nucleotideFrequencies.R -n "$sampleID" -d $output --offsetmax $value $tables
  done

done
