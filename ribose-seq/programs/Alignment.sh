#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program aligns trimmed reads to reference genome using Bowtie2 and de-duplicates reads based on UMI's
#Note: FASTQ files must be located in the users's Sequencing folder (/LocalDirectory/Ribose-Map/Sequencing)

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [-i] 'Sample(s)' [-u] 'UMI' [-m] 'Min' [-p] 'Path' [-b] 'Index' [-d] 'Directory' [-h]
		-i Input sample(s) (e.g., FS1, FS2, FS3)
		-u Length of UMI (e.g., NNNNNNNN or NNNNNNNNNNN)
		-m Minimum length of read to retain after trimming (e.g., 50)
		-p Path (e.g., /projects/home/agombolay3/data/bin/Trimmomatic-0.36)
		-b Basename of Bowtie2 index (e.g., sacCer2, pombe, ecoli, mm9, or hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "i:u:m:p:b:d:h" opt; do
    case "$opt" in
        #Allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Allow only one input argument
	u ) UMI=$OPTARG ;;
	m ) MIN=$OPTARG ;;
	p ) path=$OPTARG ;;
	b ) index=$OPTARG ;;
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
	
	#Input file
	fastq=$directory/Ribose-Map/Sequencing/$sample.fastq
	
	#Output directory and files
	output=$directory/Ribose-Map/Results/$index/$sample/Alignment
	finalBam=$output/$sample.bam; statistics=$output/Bowtie2.log
	
	#Create folder
	mkdir -p $output
	
#############################################################################################################################
	#STEP 1: Trim FASTQ files based on quality and Illumina adapter content
	java -jar $path/trimmomatic-0.36.jar SE -phred33 $fastq QCtrimmed.fastq \
	ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:$MIN
	
	#STEP 2: Extract UMI from 5' ends of reads (append UMI to read name for later)
	umi_tools extract -I QCtrimmed.fastq -p $UMI --supress-stats -S UMItrimmed.fastq
	
	#STEP 3: Reverse complement (RC) reads (R = RC of 5' base)
	cat UMItrimmed.fastq | seqtk seq -r - > reverseComplement.fastq

#############################################################################################################################
	#STEP 4: Align reads to reference and save alignment statistics file
	bowtie2 -x $index -U reverseComplement.fastq 2> $statistics > temp.sam
	
	#STEP 5: Extract mapped reads, convert SAM file to BAM, and sort and index BAM file
	samtools view -bSF4 temp.sam | samtools sort - -o temp.bam && samtools index temp.bam

#############################################################################################################################
	#STEP 6: De-duplicate reads based on UMI and start positions and sort and index BAM file
	umi_tools dedup -I temp.bam -v 0 | samtools sort - -o $finalBam && samtools index $finalBam

	#Remove temporary files
	rm -f reverseComplement.fastq temp.*
		
	#Notify user that alignment step is complete for which samples
	echo "Alignment of $sample to $index reference genome is complete"
done
