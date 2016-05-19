#!/usr/bin/env bash

#Author: Alli Gombolay
#This program automatically runs input FASTQ.GZ files through the Ribose-seq analysis pipeline
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program given as the first entry on the commnad-line (i.e., getOptions.sh))
program=$0

#Usage statement for program; will be displayed to standard output if user specifies "-h" option
function usage () {
        echo "Usage: $program [-a] 'sample1 sample2 sample3 etc.' [-b] 'basename of Bowtie index' [-h]
          -a Runs ribose-seq pipeline on input sample.fastq files using Bowtie index 
          -b Runs ribose-seq pipeline on input sample.fastq files using Bowtie index
          -h Displays help menu describing options"
}

#Use getOpts function to create command-line options (i.e., "-a", "-b", and "-h")
while getopts "a:b:h" opt; do
    case $opt in
        a ) fastq=($OPTARG) ;; #Specify input as an array to allow multiple input arguments
        b ) index=$OPTARG ;; #Specify input as a variable to allow only one input argument
        h ) usage ;; #Specify "-h" (help) option as usage statement
    esac
done

if [ "$1" = "-h" ]; then
        exit
        else 1_alignment.sh;
fi
