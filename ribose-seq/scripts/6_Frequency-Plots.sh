#! /usr/bin/env bash

#Author: Alli Gombolay
#This program plots the frequencies of the nucleotides located in the assembled samples.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program (6_Frequency-Plots.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-i] 'sample1 etc.' [-d] 'Ribose-seq directory' [-h]
        -i Sample names of input BAM files (i.e, sample1 for sample1.bam)
        -d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-d] and [-h])
while getopts "i:d:h" opt;
do
    case $opt in
    	#Specify input as arrays to allow multiple input arguments
        i ) samples=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	d ) directory=$OPTARG ;;
    	#If user specifies [-h], print usage statement
    	h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

offset_values=(100 50 15)

#modes=("all" "only-mitochondria" "no-mitochondria" "only-2micron")
modes=("no-mitochondria")

input=$directory/ribose-seq/results/$samples/nucleotideFrequencies/

output="$directory/ribose-seq/results/$samples/plots/nucleotideFrequencies"

if [[ ! -d $output ]];
then
  mkdir -p $output
fi

for index in ${!modes[@]};
do

  mode=${modes[$index]}

  for value in ${offset_values[@]};
  do
    
    sampleID="$samples.subset-$mode"
    
    tables="$input/$samples.$mode.nucleotideFrequencies.testFile5.tab"
    
    Rscript 6_Frequency-Plots.R -n "$sampleID" -d $output --offsetmax $value $tables
    
  done

done
