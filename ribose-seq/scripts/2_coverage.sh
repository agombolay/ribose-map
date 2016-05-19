#!/usr/bin/env bash

#Author: Alli Gombolay
#This script calculates the genome coverage from the alignment results.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#INPUT
chromosomeSizes=$HOME/data/ribose-seq/data/sacCer2.chromosome.sizes
finalBAM=$HOME/data/ribose-seq/results/$sample/alignment/$sample.bam

output2=$outputDirectory/$samples/bedgraphs

if [[ ! -d $output2 ]]; then
    mkdir -p $output2
fi

#OUTPUT
bgBothStrands=$HOME/data/ribose-seq/results/FS1/bedGraphs/$sample.bothStrands.coverage.bg
bgPositiveStrand=$HOME/data/ribose-seq/results/FS1/bedGraphs/$sample.positiveStrand.coverage.bg
bgNegativeStrand=$HOME/data/ribose-seq/results/FS1/bedGraphs/$sample.negativeStrand.coverage.bg

#CALCULATE GENOME COVERAGE
bedtools genomecov -ibam $finalBAM -g $chromosomeSizes -5 -bg > $bgBothStrands
bedtools genomecov -ibam $finalBAM -g $chromosomeSizes -5 -strand + -bg > $bgPositiveStrand
bedtools genomecov -ibam $finalBAM -g $chromosomeSizes -5 -strand - -bg > $bgNegativeStrand

#Explanation of options used in step above:
#"-5": Calculate coverage of only 5’ positions
#"-g": Genome file containing chromosome sizes
#"-bg": Report coverage in bedGraph file format
#"-strand": Calculate coverage of + or - strand
#"-ibam": Specify input file as BAM file format
