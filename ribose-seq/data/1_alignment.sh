#!/usr/bin/env bash

#Author: Alli Gombolay
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake
#This program removes UMI's from sequencing reads, aligns reads to reference genome, and de-duplicates reads.

for sample in ${input[@]}; do

	#VARIABLE SPECIFICATION
	#Length of UMI (Unique Molecular Identifiers)
	UMI=NNNNNNNN

	#INPUT FILES
	#Location of .fastq sequencing files
	fastq=$directory/ribose-seq/fastq/$sample.fastq

	#OUTPUT
	#Location of output directory
	output=$directory/ribose-seq/results/$sample/alignment

	#Create directory for output
	if [[ ! -d $output ]]; then
    		mkdir -p $output 
	fi

	#Location of files with trimmed UMI
	umiTrimmed=$directory/$output/$sample.umiTrimmed.fastq.gz

	#Intermediate files
	intermediateSAM=$directory/$output/$sample.intermediate.sam
	intermediateBAM=$directory/$output/$sample.intermediate.bam

	sortedBAM=$directory/$output/$sample.sorted.bam

	#Final BAM files
	finalBAM=$directory/$output/$sample.bam

	#Output file of Bowtie alignment statistics
	statistics=$directory/$output/$sample.statistics.txt

	#BED file
	BED=$directory/$output/$sample.bed.gz

	#ALIGNMENT
	
	#1. Trim UMI from input .fastq files and compress output files
	python2.7 umitools.py trim $fastq $UMI | gzip -c > $umiTrimmed

	#2. Align UMI trimmed reads to reference genome and output alignment statistics
	zcat $umiTrimmed | bowtie --uniq --sam $reference - 2> $statistics 1> $intermediateSAM

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
	samtools sort $intermeidateBAM > $sortedBAM

	#3. De-duplicate reads based on UMI's and compress BED files
	python2.7 umitools.py rmdup $sortedBAM $finalBAM | gzip -c > $BED
	
	#Remove intermediate files
	rm $intermediateSAM $intermediateBAM

done
