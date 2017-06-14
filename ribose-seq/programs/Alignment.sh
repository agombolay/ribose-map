#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program aligns trimmed reads to reference genome using Bowtie2 and de-duplicates reads based on UMI's
#Note: FASTQ files must be located in the users's FASTQ-Files folder (/LocalDirectory/Ribose-Map/FASTQ-Files)

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [options]
		-s Sample name(s) (e.g., FS1, FS2, FS3)
		-a Input Read 1 FASTQ filename (forward)
		-b Input Read 2 FASTQ filename (reverse)
		-u Length of UMI (e.g., NNNNNNNN or NNNNNNNNNNN)
		-m Minimum length of read to retain after trimming (e.g., 50)
		-p Path (e.g., /projects/home/agombolay3/data/bin/Trimmomatic-0.36)
		-t Type of Illumina Sequencing (e.g., SE = Single end, PE = Paired end)
		-i Basename of Bowtie2 index (e.g., sacCer2, pombe, ecoli, mm9, or hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:a:b:u:m:t:p:i:d:h" opt; do
    case "$opt" in
        #Allow multiple input arguments
        s ) sample=($OPTARG) ;;
	#Allow only one input argument
	a ) read1=$OPTARG ;;
	b ) read2=$OPTARG ;;
	u ) UMI=$OPTARG ;;
	m ) MIN=$OPTARG ;;
	t ) type=$OPTARG ;;
	p ) path=$OPTARG ;;
	i ) index=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#############################################################################################################################
#Align reads to reference
for sample in ${sample[@]}; do
	
	#Input files
	Read1Fastq=$directory/Ribose-Map/FASTQ-Files/$read1
	Read2Fastq=$directory/Ribose-Map/FASTQ-Files/$read2

	#Output files
	statistics=$directory/Ribose-Map/Results/$index/$sample/Alignment/Bowtie2.log
	mapped=$directory/Ribose-Map/Results/$index/$sample/Alignment/$sample-MappedReads.bam

	#Create directory
	mkdir -p $directory/Ribose-Map/Results/$index/$sample/Alignment
	
#############################################################################################################################
	#STEP 1: Trim FASTQ files based on quality and adapter content
	#Single End Reads
	if [ $type == "SE" ]; then
		java -jar $path/trimmomatic-0.36.jar SE -phred33 $Read1Fastq Read1.fq \
		ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:10 MINLEN:$MIN
	#Paired End Reads
	elif [ $type == "PE" ]; then
		java -jar $path/trimmomatic-0.36.jar PE -phred33 $Read1Fastq $Read2Fastq Read1.fq Unpaired1.fq \
		Read2.fq Unpaired2.fq ILLUMINACLIP:$path/adapters/TruSeq3-PE.fa:2:30:10 TRAILING:10 MINLEN:$MIN
	fi
	
	#STEP 3: Reverse complement reads (Ribo = RC of 5' base of read)
	#Single End Reads
	if [ $type == "SE" ]; then
		cat Read1.fq | seqtk seq -r - > R1Reverse.fq
	#Paired End Reads
	elif [ $type == "PE" ]; then
		cat Read1.fq | seqtk seq -r - > R1Reverse.fq
		cat Read2.fq | seqtk seq -r - > R2Reverse.fq
	fi
	
	#STEP 2: Extract UMI from 5' ends of reads (append UMI to read name)
	umi_tools extract -I R1Reverse.fq -p $UMI --3prime --supress-stats -S R1Trimmed.fq
	
#############################################################################################################################
	#STEP 4: Align reads to reference genome and save Bowtie2 statistics to file
	#Single End Reads
	if [ $type == "SE" ]; then
		bowtie2 -x $index -U R1Trimmed.fq 2> $statistics > temp.sam
	#Paired End Reads
	elif [ $type == "PE" ]; then
		bowtie2 -x $index -1 R1Trimmed.fq -2 R2Reverse.fq 2> $statistics -S temp.sam
	fi
	
	#STEP 5: Extract mapped reads, convert SAM file to BAM, and sort/index BAM file
	#Single End Reads
	if [ $type == "SE" ]; then
		samtools view -bSF4 temp.sam | samtools sort - -o temp.bam; samtools index temp.bam
	#Paired End Reads
	elif [ $type == "PE" ]; then
		samtools view -bSf66 temp.sam | samtools sort - -o temp.bam; samtools index temp.bam
	fi
	
#############################################################################################################################
	#STEP 6: Remove PCR duplicates based on UMI and position and sort/index BAM file
	umi_tools dedup -I temp.bam -v 0 | samtools sort - -o $mapped; samtools index $mapped

	#Notify user that alignment step is complete for which samples
	echo "Alignment of $sample to $index reference genome is complete"
	
	#Remove temporary files
	rm -f QCtrimmed.fastq UMItrimmed.fastq reverseComplement.fastq temp.*

done
