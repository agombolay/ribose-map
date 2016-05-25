#! /usr/bin/env bash

#Author: Alli Gombolay
#This script calculates nucleotide frequencies in sacCer2 genome (chr I-XVI, chr M, and 2micron).
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#INPUT
fasta=$HOME/data/ribose-seq/data/reference/sacCer2.fa

#OUTPUT

output="$directory/ribose-seq/results/backgroundNucleotideFrequencies"

if [[ ! -d $output ]];
then
  mkdir $output
fi

#CALCULATION OF NUCLEOTIDE FREQUENCIES

#Calculate frequencies of nucleotides in sacCer2 genome
output="$HOME/data/ribose-seq/results/nucleotideFrequencies/genome.nuc.freqs.tab"
python2.7 4_nucleotideFrequencies.py $fasta --region-size-minimum 1 --region-size-maximum 3 --verbose > $output

#Calculate frequencies of nucleotides in chrM
output="$HOME/data/ribose-seq/results/nucleotideFrequencies/chrM.nuc.freqs.tab"
python2.7 4_nucleotideFrequencies.py $fasta --region-size-minimum 1 --region-size-maximum 3 --only-chrom chrM --verbose > $output

#Calculate frequencies of nucleotides in 2micron 
output="$HOME/data/ribose-seq/results/nucleotideFrequencies/2micron.nuc.freqs.tab"
python2.7 4_nucleotideFrequencies.py $fasta --region-size-minimum 1 --region-size-maximum 3 --only-chrom 2micron --verbose > $output
