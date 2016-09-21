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

for sample in ${fastq[@]}; do
	#Extract names from filepaths
	filename=$(basename "${sample}")
	sample="${filename%.*}"
	
	#Extract directory from filepaths
	directory0=$(dirname "${fastq}")
	
	#Location of FASTQ files
	fastq=$directory0/$sample.fastq

	echo $fastq
done

#Align FASTQ files to reference genome
for sample in ${fastq[@]}; do

	#Extract names from filepaths
	filename=$(basename "${sample}")
	sample="${filename%.*}"
	
	#Extract directory from filepaths
	directory0=$(dirname "${fastq}")
	
	#Location of FASTQ files
	fastq=$directory0/$sample.fastq
	
	#OUTPUT
	#Location of output directory
	output=$directory/ribose-seq/results/$index/$sample/Alignment/

	#Create directory if not present
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

	#Length of Unique Molecular Identifiers (UMI)
	UMI=NNNNNNNN

#############################################################################################################################
	#1. Reverse complement reads
	seqtk seq -r $fastq > $reverseComplement

	#2. Trim UMI from 3' ends of reads; compress file
	umitools trim --end 3 $reverseComplement $UMI | gzip -c > $umiTrimmed

	#3. Align reads to reference genome with Bowtie version 1 or 2
	#Bash: "-": standard input; "2>": Redirect standard error; "1>": Redirect standard output
	#Bowtie: "-m 1": Return only unique reads; "--sam": Print alignment results in SAM format
	if [ $version == "1" ]; then
		zcat $umiTrimmed | bowtie -m 1 --sam $index - 2> $statistics 1> $intermediateSAM
	elif [ $version == "2" ]; then
		zcat $umiTrimmed | bowtie2 -x $index - 2> $statistics 1> $intermediateSAM
	fi

	#4. Convert SAM file to BAM
	#SAMtools: #"-S": Input format is SAM; "-h": Include header in output;
	#"-u": Output as uncompressed BAM; #"-F4": Do not output unmapped reads
	samtools view -ShuF4 $intermediateSAM > $intermediateBAM
	
	#5. Sort intermediate BAM files
	samtools sort $intermediateBAM > $sortedBAM

	#6. Index sorted BAM files
	samtools index $sortedBAM
	
	#7. De-duplicate reads based on UMIs; compress file
	umitools rmdup $sortedBAM $finalBAM | gzip -c > $BED
	
	#8. Index final BAM files
	samtools index $finalBAM
	
done

echo "Alignment of $sample to reference genome is complete"
