#!/usr/bin/env bash

#Author: Alli Gombolay
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake
#This program removes UMI's from sequencing reads, aligns reads to reference genome, and deduplicates reads.

#SPECIFY COMMAND LINE OPTIONS

program=$0

function usage () {
        echo "Usage: $program -g 'input1 input2 input3 etc.' -r 'basename of Bowtie index'"
}

while getopts "g:r:h" opt; do
    case $opt in
        g ) input=($OPTARG) #Specify input as an array
            echo "Aligning ${input[@]} genomes to reference" ;;
	r ) reference=($OPTARG) ;; #Specify input as an array
        h ) usage ;;
    esac
done

for sample in ${input[@]}; do

#VARIABLE SPECIFICATION

#Path to umitools
umitoolsPath=/projects/home/agombolay3/.local/lib/python2.7/site-packages/umitools-2.1.1-py2.7.egg/umitools/

#Length of UMI's
UMI=NNNNNNNN

#Align modes
ALIGN_MODES=("all")

#Arguments of alignment
ALIGN_ARGUMENTS=("--all")

#INPUT

#Location of raw FASTQ files
unprocessedFASTQ=$HOME/data/sequencingResults/run2_February2016/compressedFiles/$sample.trimmed.fastq.gz

#OUTPUT

#LOCATION OF OUTPUT FILES
output=$HOME/data/ribose-seq/results/$sample/alignment

#CREATE DIRECTORY STRUCTURE FOR OUTPUT FILES
if [[ ! -d $output ]]; then
    mkdir -p $output 
fi

#Location of files with trimmed UMI
umiTrimmed=$output/$sample.umiTrimmed.fastq.gz

#Intermediate files
SAM=$output/$sample.sam
BAM=$output/$sample.bam

sortedBAM=$output/$sample.sorted.bam

#Main output BAM files
finalBAM=$output/$sample.final.bam

#Output file detailing Bowtie statistics
statistics=$output/$sample.statistics.txt

#Final output BED file
finalBED=$output/$sample.bed.gz

#ALIGNMENT
#1. Trim UMI from FASTQ files and then compress output files
python2.7 $umitoolsPath/umitools.py trim $unprocessedFASTQ $UMI | gzip -c > $umiTrimmed

#2. Align FASTQ reads with Bowtie and output alignment statistics
zcat $umiTrimmed | bowtie --all --sam $reference - 2> $statistics 1> $SAM

#Explanation of options used in step above:
#"-": standard input
#2>: Redirect standard error to file
#1>: Redirect standard output to file
#"--all": Return all valid alignments per read
#"--sam": Print alignment results in SAM format
#reference = Basename of Bowtie index to be searched

#Convert SAM file to BAM file
samtools view -bS $SAM > $BAM

#Explanation of options used in step above:
#"-b": Output in BAM format
#"-S": Input in SAM format

#Sort BAM file
samtools sort $BAM > $sortedBAM

#3. De-duplicate reads based on UMI information
python2.7 $umitoolsPath/umitools.py rmdup $sortedBAM $finalBAM | gzip -c > $finalBED;

done
