#!/usr/bin/env bash

#Author: Alli Gombolay
#This program removes UMI's from reads, aligns reads to reference genome, and de-duplicates reads
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program (1_Alignment.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-i] '/path/to/file1.fastq etc.' [-b] 'Bowtie index' [-d] 'Ribose-seq directory' [-h]
          -i Filepaths of input FASTQ files 
          -b Basename of Bowtie index to be searched
          -d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-i], [-b], [-d], and [-h])
while getopts "i:b:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) fastq=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	b ) index=$OPTARG ;;
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

#Align all of the input sequence data to the specified reference genome
for samples in ${fastq[@]};
do
	
	#Extract sample names from filepaths
	filename=$(basename "${samples}")
	samples="${filename%.*}"
	
	#Extract input directories from filepaths
	inputDirectory=$(dirname "${fastq}")
	
	#VARIABLE SPECIFICATION
	#Length of UMI (Unique Molecular Identifiers)
	UMI=NNNNNNNN

	#INPUT
	#Location of raw sequencing files
	input=$inputDirectory/$samples.fastq

	#OUTPUT
	#Location of output "ribose-seq" alignment directory
	output=$directory/ribose-seq/results/$index/$samples/Alignment/

	#Create directory for output if it does not already exist
	if [[ ! -d $output ]];
	then
    		mkdir -p $output
	fi

	#Location of reverse complement files of input FASTQ files
	reverseComplement=$output/$samples.reverse.complement.fastq

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

	#Obtain reverse complement of input FASTQ files
	seqtk seq -r $input > $reverseComplement
	
	#1. Trim UMI from 3' ends of reads and compress output files
	umitools trim --end 3 $reverseComplement $UMI | gzip -c > $umiTrimmed

	#2. Align UMI trimmed reads to reference genome and output alignment statistics
	#zcat $umiTrimmed | bowtie -m 1 --sam $index --un $samples.unmappedReads --max $samples.extraReads - 2> $statistics 1> $intermediateSAM
	
	#zcat $umiTrimmed | bowtie -m 1 --sam $index - 2> $statistics 1> $intermediateSAM
	zcat $umiTrimmed | bowtie2 -x $index - 2> $statistics 1> $intermediateSAM

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
	
	#3. De-duplicate reads based on UMI's and compress BED files
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
