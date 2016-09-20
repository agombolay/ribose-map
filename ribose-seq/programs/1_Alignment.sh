#!/usr/bin/env bash

#Author: Alli Gombolay
#This program removes UMI's from reads, aligns reads to reference genome, and de-duplicates reads
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 1_Alignment.sh [-i] 'FASTQ' [-b] 'Index' [-v] 'Bowtie Version' [-d] 'Directory' [-h]
		-i Filepaths of input FASTQ files ('/path/to/file1.fastq' etc.) 
		-b Basename of Bowtie index to be searched (sacCer2, chrM, ecoli, hg38, etc.)
		-v Version of Bowtie program to use (Version 1 = Bowtie1; Version 2 = Bowtie2)
		-d Local directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Use getopts function to create the command-line options ([-i], [-b], [-d], [-v], and [-h])
while getopts "i:b:d:v:h" opt; do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) fastq=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	b ) index=$OPTARG ;;
	v ) version=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Align FASTQ files to reference genome
for sample in ${fastq[@]}; do
	
	#Extract sample names from filepaths
	filename=$(basename "${sample}")
	sample="${filename%.*}"
	
	#Extract input directory from filepaths
	inputDirectory=$(dirname "${fastq}")
	
	#VARIABLE SPECIFICATION
	#Length of Unique Molecular Identifiers (UMI)
	UMI=NNNNNNNN

	#INPUT
	#Location of FASTQ files
	input=$inputDirectory/$sample.fastq

	#OUTPUT
	#Location of output directory
	output=$directory/ribose-seq/results/$index/$sample/Alignment/

	#Create directory if it does not already exist
	if [[ ! -d $output ]]; then
    		mkdir -p $output
	fi

	#Location of reverse complemented FASTQ files 
	reverseComplement=$output/$sample.reverse-complement.fastq

	#Location of UMI trimmed FASTQ files
	umiTrimmed=$output/$sample.UMI-trimmed.fastq.gz

	#Intermediate files
	intermediateSAM=$output/$sample.intermediate.sam
	intermediateBAM=$output/$sample.intermediate.bam
	sortedBAM=$output/$sample.sorted.bam

	#Final BAM files
	finalBAM=$output/$sample.bam

	#File containing Bowtie alignment statistics
	statistics=$output/$sample.Alignment-Statistics.txt

	#BED file
	BED=$output/$samples.bed.gz

	#ALIGNMENT

	#Reverse complement reads
	seqtk seq -r $input > $reverseComplement
	
	#1. Trim UMI from 3' ends of reads; compress file
	umitools trim --end 3 $reverseComplement $UMI | gzip -c > $umiTrimmed

	#2. Align UMI trimmed reads to reference genome and output alignment statistics
	#zcat $umiTrimmed | bowtie -m 1 --sam $index --un $samples.unmappedReads --max $samples.extraReads - 2> $statistics 1> $intermediateSAM
	
	#Align reads with Bowtie version 1 or Bowtie version 2
	if [ $version == "1" ]; then
		zcat $umiTrimmed | bowtie -m 1 --sam $index - 2> $statistics 1> $intermediateSAM
	
	elif [ $version == "2" ]; then
		zcat $umiTrimmed | bowtie2 -x $index - 2> $statistics 1> $intermediateSAM
	fi

	#Bash functions used above:
	#"-": standard input
	#2>: Redirect standard error to file
	#1>: Redirect standard output to file
	
	#Bowtie options used above:
	#"-m 1": Return only unique reads
	#"--sam": Print alignment results in SAM format
	#reference: Basename of Bowtie index to be searched
	
	#Convert SAM files to BAM files
	#samtools view -bS $intermediateSAM > $intermediateBAM
	samtools view -ShuF4 $intermediateSAM > $intermediateBAM
	
	#SAMtools options used above:
	#"-S": Input in SAM format
	#"-h": Include header in output
	#"-u": Output as uncompressed BAM
	#"-F4": Do not output unmapped reads

	#Sort intermediate BAM files
	samtools sort $intermediateBAM > $sortedBAM

	#Index sorted BAM files
	samtools index $sortedBAM
	
	#3. De-duplicate reads based on UMIs; compress file
	umitools rmdup $sortedBAM $finalBAM | gzip -c > $BED
	
	#Index final BAM files
	samtools index $finalBAM
	
	#Ask the user if he/she would like to delete all of the intermediate files
	#echo "Would you like to delete all of the intermediate files for $samples? [y/n]"
	
	#"answer" is the variable representing the user's answer
	#read answer
	
	#If the user answers "y," then remove the specified files
	#if [ $answer == "y" ];
        #then
        #	rm $umiTrimmed $intermediateSAM $intermediateBAM $sortedBAM;
        #fi

done
