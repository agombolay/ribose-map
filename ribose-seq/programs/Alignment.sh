#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program aligns trimmed reads to reference genome using Bowtie2 and de-duplicates reads based on UMI's
#Note: FASTQ files must be located in the users's FASTQ-Files folder (/LocalDirectory/Ribose-Map/FASTQ-Files)

#Usage statement
function usage () {
echo "Usage: Alignment.sh [-s] 'Sample(s)' [-u] 'UMI' [-m] 'Min' [-p] 'Path' [-t] 'Type' [-i] 'Index' [-d] 'Directory' [-h]
	-s Sample name(s) (e.g., FS1, FS2, FS3)
	-a Input Read 1 FASTQ.GZ filename (forward)
	-b Input Read 2 FASTQ.GZ filename (reverse)
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
	a ) fastq1=$OPTARG ;;
	b ) fastq2=$OPTARG ;;
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
	fastq1=$directory/Variant-Calling/Sequencing/$Read1
	fastq2=$directory/Variant-Calling/Sequencing/$Read2

	#Output files
	statistics=$directory/Ribose-Map/Results/$index/$sample/Alignment/Bowtie2.log
	mapped=$directory/Ribose-Map/Results/$index/$sample/Alignment/$sample-MappedReads.bam

	#Create directory
	mkdir -p $directory/Ribose-Map/Results/$index/$sample/Alignment
	
#############################################################################################################################
	#STEP 1: Trim FASTQ files based on quality and Illumina adapter content
	#Single End Reads
	if [ $type == "SE" ]; then
		java -jar $path/trimmomatic-0.36.jar SE -phred33 $fastq1 QCtrimmed.fastq \
		ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:10 MINLEN:$MIN
	#Paired End Reads
	elif [ $type == "PE" ]; then
		java -jar $path/trimmomatic-0.36.jar PE -phred33 $fastq1 $fastq2 R1Paired.fq.gz R1Unpaired.fq.gz \
		R2Paired.fq.gz R2Unpaired.fq.gz ILLUMINACLIP:$path/adapters/TruSeq3-PE.fa:2:30:10 TRAILING:10 MINLEN:$MIN
	fi
	
	#STEP 2: Extract UMI from 5' ends of reads (append UMI to read name for later)
	umi_tools extract -I QCtrimmed.fastq -p $UMI --supress-stats -S UMItrimmed.fastq
	
	#STEP 3: Reverse complement (RC) reads (R = RC of 5' base)
	cat UMItrimmed.fastq | seqtk seq -r - > reverseComplement.fastq
	
#############################################################################################################################
	#STEP 4: Align reads to reference genome and save Bowtie2 log file
	bowtie2 -x $index -U reverseComplement.fastq 2> $statistics > temp.sam
	
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
