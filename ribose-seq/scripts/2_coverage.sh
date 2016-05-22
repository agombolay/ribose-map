#!/usr/bin/env bash

#Author: Alli Gombolay
#This script calculates the genome coverage based on the alignment results from 1_alignment.sh
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program (1_alignment.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-i] 'sample1 sample2 sample3 etc.' [-r] 'reference' [-d] 'directory' [-h]
          -i Sample names of input BAM files (i.e, FS1 for FS1.bam)
          -r File containing sizes in base pairs of chromosomes
          -d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-a], [-b], [-o], and [-h])
while getopts "i:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) list=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	r ) reference=$OPTARG ;;
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

for samples in ${fastq[@]};
do

	#INPUT
	#Location of BAM files
	input=$directory/ribose-seq/results/$samples/alignment/$samples.bam
	
	#Location of files containing sizes in base pairs of chromosomes
	chromosomeSizes=$directory/ribose-seq/data/$reference.chromosome.sizes
	
	#OUTPUT
	output=$directory/ribose-seq/results/$samples/bedgraphs

	if [[ ! -d $output ]];
	then
		mkdir -p $output
	fi

	#Location of output bedgraph files containing genome coverage information
	BothStrands=$directory/ribose-seq/results/$samples/bedGraphs/$samples.bothStrands.coverage.bg
	PositiveStrands=$directory/ribose-seq/results/$samples/bedGraphs/$samples.positiveStrands.coverage.bg
	NegativeStrands=$directory/ribose-seq/results/$samples/bedGraphs/$samples.negativeStrands.coverage.bg

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
