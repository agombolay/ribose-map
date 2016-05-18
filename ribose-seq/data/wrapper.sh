#!/usr/bin/env bash

#Author: Alli Gombolay
#This program automatically runs input FASTQ.GZ files through the Ribose-seq analysis pipeline
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#SPECIFY COMMAND LINE OPTIONS

program=$0

function usage () {
        echo "Usage: $program [-a] 'input1 input2 input3 etc.' [-b] 'basename of Bowtie index' [-h]
          -a Runs ribose-seq pipeline on input .fastq files using Bowtie index 
          -b Runs ribose-seq pipeline on input .fastq files using Bowtie index
          -h Displays help menu describing options"
}

while getopts "a:b:h::" opt; do
    case $opt in
        g ) input=($OPTARG) #Specify input as an array
            echo "Aligning ${input[@]} genomes to reference" ;;
        r ) reference=($OPTARG) ;; #Specify input as an array
        h ) usage ;;
    esac
done
