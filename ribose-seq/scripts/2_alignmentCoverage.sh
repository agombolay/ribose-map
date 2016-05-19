#!/usr/bin/env bash

#Author: Alli Gombolay
#This script calculates the genome coverage from the alignment results.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#INPUT
#chromosomeSizes=$HOME/data/ribose-seq/data/reference/sacCer2.chrom.sizes
chromosomeSizes=$HOME/data/ribose-seq/data/sacCer2.chrom.sizes
finalBAM=$HOME/data/ribose-seq/results/$sample/alignment/$sample.final.bam

output2=$outputDirectory/$samples/bedgraphs

if [[ ! -d $output2 ]]; then
    mkdir -p $output2
fi

#OUTPUT BEDGRAPH FILES
bedGraphBothStrands=$HOME/data/ribose-seq/results/FS1/bedGraphs/$sample.bothStrands.coverage.bg
bedGraphPositiveStrand=$HOME/data/ribose-seq/results/FS1/bedGraphs/$sample.positiveStrand.coverage.bg
bedGraphNegativeStrand=$HOME/data/ribose-seq/results/FS1/bedGraphs/$sample.negativeStrand.coverage.bg

#CALCULATE GENOME COVERAGE
bedtools genomecov -5 -bg -g $CHROM_SIZES -ibam $finalBAM > $bedGraphBothStrands
bedtools genomecov -5 -strand + -bg -g $CHROM_SIZES -ibam $finalBAM  > $bedGraphPositiveStrand
bedtools genomecov -5 -strand - -bg -g $CHROM_SIZES -ibam $finalBAM > $bedGraphNegativeStrand

#Explanation of options used in step above:
#"-5": Calculate coverage of only 5’ positions
#"-bg": Report coverage in bedGraph file format
#"-strand": Calculate coverage of + or - strand
#"-ibam": Specify input file as BAM file format
