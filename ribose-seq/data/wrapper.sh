#!/usr/bin/env bash

#Author: Alli Gombolay
#This program automatically runs input FASTQ.GZ files through the Ribose-seq analysis pipeline
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program given as the first entry on the commnad-line (i.e., getOptions.sh))
program=$0

#Usage statement for program; will be displayed to standard output if user specifies "-h" option
function usage () {
        echo "Usage: $program [-a] 'input1 input2 input3 etc.' [-b] 'basename of Bowtie index' [-h]
          -a Runs ribose-seq pipeline on input .fastq files using Bowtie index 
          -b Runs ribose-seq pipeline on input .fastq files using Bowtie index
          -h Displays help menu describing options"
}

#Use getOpts function to create command-line options (i.e., "-a", "-b", and "-h")
while getopts "a:b:h::" opt; do
    case $opt in
        g ) fastq=($OPTARG) ;; #Specify input as an array to allow multiple input arguments
        r ) index=($OPTARG) ;; #Specify input as an array to allow multiple input arguments
        h ) usage ;; #Specify "-h" (help) option as usage statement
    esac
done
