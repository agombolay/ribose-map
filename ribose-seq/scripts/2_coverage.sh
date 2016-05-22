#!/usr/bin/env bash

#Author: Alli Gombolay
#This script calculates the genome coverage based on the alignment results from 1_alignment.sh
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

for samples in ${fastq[@]};
do

    #INPUT
	#Location of BAM files
	input=$outputDirectory/ribose-seq/results/$samples/alignment/$samples.bam
	
	#Location of files containing sizes in base pairs of chromosomes
	chromosomeSizes=$outputDirectory/ribose-seq/data/$reference.chromosome.sizes
	
	#OUTPUT
    output=$outputDirectory/ribose-seq/results/$samples/bedgraphs

    if [[ ! -d $output ]];
    then
        mkdir -p $output
    fi

    #Location of output bedgraph files
    BothStrands=$outputDirectory/ribose-seq/results/$samples/bedGraphs/$samples.bothStrands.coverage.bg
    PositiveStrands=$outputDirectory/ribose-seq/results/$samples/bedGraphs/$samples.positiveStrands.coverage.bg
    NegativeStrands=$outputDirectory/ribose-seq/results/$samples/bedGraphs/$samples.negativeStrands.coverage.bg

    #CALCULATE GENOME COVERAGE
    bedtools genomecov -ibam $input -g $chromosomeSizes -5 -bg > $BothStrands
    bedtools genomecov -ibam $input -g $chromosomeSizes -5 -strand + -bg > $PositiveStrands
    bedtools genomecov -ibam $input -g $chromosomeSizes -5 -strand - -bg > $NegativeStrands

    #bedtools options used above:
    #"-5": Calculate coverage of only 5’ positions
    #"-g": Genome file containing chromosome sizes
    #"-bg": Report coverage in bedGraph file format
    #"-strand": Calculate coverage of + or - strand
    #"-ibam": Specify input file as BAM file format

done
