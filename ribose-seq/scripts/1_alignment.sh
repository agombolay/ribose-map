#!/usr/bin/env bash

#Author: Alli Gombolay
#This program removes UMI's from reads, aligns reads to reference genome, and de-duplicates reads
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program given as the first entry on the commnad-line (i.e., getOptions.sh))
program=$0

#Usage statement for program; will be displayed to standard output if user specifies "-h" option
function usage () {
        echo "Usage: $program [-a] 'filepath1 filepath2 etc.' [-b] 'basename of Bowtie index' [-o] 'outputDirectory' [-h]
          -a Runs ribose-seq pipeline on input sample.fastq files using Bowtie index 
          -b Runs ribose-seq pipeline on input sample.fastq files using Bowtie index
          -o Saves all results to the specified output directory
          -h Displays help menu describing options"
}

#Use getOpts function to create command-line options (i.e., "-a", "-b", "-o," and "-h")
while getopts "a:b:o:h" opt; do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        a ) fastq=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	b ) index=$OPTARG ;;
	o ) outputDirectory=$OPTARG ;;
        #Specify "-h" (help) option as usage statement written above
        h ) usage ;;
    esac
done

path=/projects/home/agombolay3/.local/lib/python2.7/site-packages/umitools-2.1.1-py2.7.egg/umitools/

for samples in ${fastq[@]}; do
	
	#Extract sample names from filepaths
	filename=$(basename "$fastq")
	samples="${filename%.*}"
	
	#Extract input directories from filepaths
	inputDirectory=$(dirname "${fastq}")
	
	#VARIABLE SPECIFICATION
	#Length of UMI (Unique Molecular Identifiers)
	UMI=NNNNNNNN

	#INPUT FILES
	#Location of raw sequencing files
	input=$inputDirectory/$samples.fastq

	#OUTPUT
	#Location of output "ribose-seq" alignment directory
	output=$outputDirectory/ribose-seq/results/$samples/alignment

	#Create directory for output
	if [[ ! -d $output ]]; then
    		mkdir -p $output
	fi

	#Location of files with trimmed UMI
	umiTrimmed=$output/$samples.umiTrimmed.fastq.gz

	#Intermediate files
	intermediateSAM=$output/$samples.intermediate.sam
	intermediateBAM=$output/$samples.intermediate.bam

	sortedBAM=$output/$samples.sorted.bam

	#Final BAM files
	finalBAM=$output/$samples.bam

	#Output file of Bowtie alignment statistics
	statistics=$output/$samples.statistics.txt

	#BED file
	BED=$output/$samples.bed.gz

	#ALIGNMENT
	
	#1. Trim UMI from raw reads and compress output files
	python2.7 $path/umitools.py trim $input $UMI | gzip -c > $umiTrimmed

	#2. Align UMI trimmed reads to reference genome and output alignment statistics
	zcat $umiTrimmed | bowtie -m 1 --sam $index - 2> $statistics 1> $intermediateSAM

	#Bash functions used above:
	#"-": standard input
	#2>: Redirect standard error to file
	#1>: Redirect standard output to file
	
	#Bowtie options used above:
	#"-m 1": Return only unique reads
	#"--sam": Print alignment results in SAM format
	#reference: Basename of Bowtie index to be searched
	
	#Convert SAM files to BAM files
	samtools view -bS $intermediateSAM > $intermediateBAM

	#SAMtools options used above:
	#"-b": Output in BAM format
	#"-S": Input in SAM format

	#Sort intermediate BAM files
	samtools sort $intermediateBAM > $sortedBAM

	#3. De-duplicate reads based on UMI's and compress BED files
	python2.7 $path/umitools.py rmdup $sortedBAM $finalBAM | gzip -c > $BED
	
	#Ask the user if he/she would like to delete all of the intermediate files
	echo "Would you like to delete all of the intermediate files for $samples? [y/n]"
	
	#"answer" is the variable representing the user's answer
	read answer
	
	#If the user answers "y," then remove the specified files
	if [ $answer == "y" ];
        	then rm $umiTrimmed $intermediateSAM $intermediateBAM $sortedBAM;
        fi

done
