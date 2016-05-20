#!/usr/bin/env bash

#Author: Alli Gombolay
#This program removes UMI's from reads, aligns reads to reference genome, and de-duplicates reads
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

path=/projects/home/agombolay3/.local/lib/python2.7/site-packages/umitools-2.1.1-py2.7.egg/umitools/

for samples in ${fastq[@]}; do

	filename=$(basename "$fastq")
	extension="${filename##*.}"
	samples="${filename%.*}"
	
	inputDirectory=$(dirname "${fastq}")
	outputDirectory=/projects/home/agombolay3/data

	#VARIABLE SPECIFICATION
	#Length of UMI (Unique Molecular Identifiers)
	UMI=NNNNNNNN

	#INPUT FILES
	#Location of raw sequencing files
	input=$inputDirectory/$samples.fastq

	#OUTPUT
	#Location of output "ribose-seq" alignment directory
	output1=$outputDirectory/ribose-seq/results/$samples/alignment

	#Create directory for output
	if [[ ! -d $output1 ]]; then
    		mkdir -p $output1 
	fi

	#Location of files with trimmed UMI
	umiTrimmed=$output1/$samples.umiTrimmed.fastq.gz

	#Intermediate files
	intermediateSAM=$output1/$samples.intermediate.sam
	intermediateBAM=$output1/$samples.intermediate.bam

	sortedBAM=$output1/$samples.sorted.bam

	#Final BAM files
	finalBAM=$output1/$samples.bam

	#Output file of Bowtie alignment statistics
	statistics=$output1/$samples.statistics.txt

	#BED file
	BED=$output1/$samples.bed.gz

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
