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
		-d Directory (e.g., /projects/home/agombolay3/data/repository/Ribose-seq-Project)"
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
	
	#Input FASTQ files
	fastq=$directory/Sequencing-Results/$sample.fastq
	
	#Output directory and create directory
	output=$directory/ribose-seq/results/$index/$sample/Alignment; mkdir -p $output

	#Output files
	finalBAM=$output/$sample-mappedReads.bam; statistics=$output/$sample-Statistics.txt
	
#############################################################################################################################
	#STEP 1: QUALITY TRIMMING
	
	#Trim FASTQ files based on quality and Illumina adapter content
	#java -jar $path/trimmomatic-0.36.jar SE -phred33 $fastq $output/QCtrimmed.fastq \
	#ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:10 SLIDINGWINDOW:5:15 MINLEN:$MIN
	
	#STEP 2: EXTRACT UMI FROM READS
	#Trim UMI from 5' ends of reads (append UMI to read name for subsequent de-duplication step)
	#umi_tools extract -I $output/QCtrimmed.fastq -p $UMI --supress-stats -S $output/UMItrimmed.fastq
	
	#STEP 3: REVERSE COMPLEMENT READS
	#Reverse complement reads (rNMP=reverse complement of 5' base)
	#cat $output/UMItrimmed.fastq | seqtk seq -r - > $output/reverseComplement.fastq
	
	#STEP 4: ALIGN READS TO REFERENCE GENOME
	#Align reads to reference using Bowtie2 and output alignment statistics
	#bowtie2 -x $index -U $output/reverseComplement.fastq -S $output/temp.sam 2> $statistics
	
	bowtie2 -x $index -U $output/reverseComplement.fastq 2> $statistics | samtools view -bSF4 - | \
	samtools sort - -o $output/mapped.bam && samtools index $output/mapped.bam
	
	#STEP 5: CONVERT SAM FILE TO BAM FILE AND SORT/INDEX IT
	#Directly convert SAM file to sorted BAM file (Save only mapped reads) and create index for BAM file
	#samtools view -bSF4 $output/temp.sam | samtools sort - -o $output/mapped.bam && samtools index $output/mapped.bam
	
	#STEP 6: DE-DUPLICATE READS BASED ON UMI AND SORT/INDEX BAM FILE
	#De-duplicate reads based on UMI and coordinates, sort BAM file, and create an index file for it
	umi_tools dedup -I $output/mapped.bam -v 0 | samtools sort - -o $finalBAM && samtools index $finalBAM

	#Remove temporary files
	rm -f $output/reverseComplement.fastq $output/temp.bam \
	$output/temp.bam.bai $output/temp.sam $output/mapped.bam
		
	#Notify user that alignment step is complete for which samples
	echo "Alignment of $sample to $index reference genome is complete"
done
