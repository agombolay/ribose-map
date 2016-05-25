#! /usr/bin/env bash

#Author: Alli Gombolay
#This script calculates nucleotide frequencies in sacCer2 genome (chr I-XVI, chr M, and 2micron).
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#INPUT
#Location of FASTA file containing sequences of sacCer2 genome
fasta=$directory/reference/sacCer2.fa

#OUTPUT

#Location of output "ribose-seq" backgroundNucleotideFrequencies directory
output="$directory/ribose-seq/results/backgroundNucleotideFrequencies"

#Create output directory if it does not already exist
if [[ ! -d $output ]];
then
  mkdir $output
fi

#CALCULATION OF BACKGROUND NUCLEOTIDE FREQUENCIES

#Calculate frequencies of nucleotides in sacCer2 genome
output1="$output/genome.nucleotide.frequencies.tab"
python2.7 4_nucleotideFrequencies.py $fasta --region-size-minimum 1 --region-size-maximum 3 --verbose > $output1

#Calculate frequencies of nucleotides in chrM
output2="$output/chrM.nucleotide.frequencies.tab"
python2.7 4_nucleotideFrequencies.py $fasta --region-size-minimum 1 --region-size-maximum 3 --only-chrom chrM --verbose > $output2

#Calculate frequencies of nucleotides in 2micron 
output3="$output/2micron.nucleotide.frequencies.tab"
python2.7 4_nucleotideFrequencies.py $fasta --region-size-minimum 1 --region-size-maximum 3 --only-chrom 2micron --verbose > $output3
