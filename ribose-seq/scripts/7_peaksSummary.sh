#! /usr/bin/env bash

#Author: Alli Gombolay
#This script generates summary tables of the peaks calling resulting from the "MACS2" program
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

output=$directory/ribose-seq/results/peaksSummary

if [[ ! -d $output ]];
then
    mkdir $output
fi

strands=("positive" "negative")

for strand in ${strands[@]};
do

        samples=""
        names=""
        
        for sample in ${SAMPLES[@]}; do
            input=$directory/ribose-seq/results/$samples/peaks
            
            experiment="$sample.strand.$strand"
            
            peakbed="$input/$experiment""_peaks.narrowPeak"
            
            samples="$samples $peakbed"
            
            names="$names $sample"
        done

        tables="$output/peaksSummary.$strand.strand.tab"

        bedtools multiinter -i $samples -names $names -g $chromosomeSizes > $tables
        
    done
done
