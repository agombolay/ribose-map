#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program raligns trimmed reads to reference genome using Bowtie2 and de-duplicates reads based on UMI's
#Note: FASTQ files must be located in the users's Sequencing-Results folder (/LocalDirectory/Sequencing-Results)

#Usage statement
function usage () {
	echo "Usage: 1_Alignment.sh [-i] 'FASTQ' [-b] 'Index' [-d] 'Directory' [-h]
		-i Sample name(s) (FS1, FS2, FS3 etc.) 
		-b Basename of Bowtie2 index to be searched (sacCer2, ecoli, mm9, hg38, etc.)
		-d Local directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Command-line options
while getopts "i:b:d:h" opt; do
    case "$opt" in
        #Allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Allow only one input argument
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
	
	#Output directory
	output=$directory/ribose-seq/results/$index/$sample/Alignment/

	#Create directory
    	mkdir -p $output

	#Intermediate files
	umiTrimmed=$output/$sample.UMI-trimmed.fastq.gz; intermediateSAM=$output/$sample.intermediate.sam;
	sortedBAM=$output/$sample.sorted.bam; reverseComplement=$output/$sample.reverse-complement.fastq

	#UMI length
	UMI=NNNNNNNNNNN
	
	#BED file
	BED=$output/$sample.bed.gz
	
	#Final BAM file
	finalBAM=$output/$sample.bam

	#File containing Bowtie2 alignment statistics
	statistics=$output/$sample.Alignment-Statistics.txt

#############################################################################################################################
	#Reverse complement reads
	seqtk seq -r $fastq > $reverseComplement

	#Trim UMI from 3' ends of reads (add UMI into read name)
	umitools trim --end 3 $reverseComplement $UMI | gzip -c > $umiTrimmed

	#Align reads to reference genome using Bowtie
	#zcat $umiTrimmed | bowtie -m 1 $index - -S $intermediateSAM 2> $statistics
	zcat $umiTrimmed | bowtie2 -x $index -U - -S $intermediateSAM 2> $statistics
	
	#Convert SAM file to BAM and sort intermediate BAM file
	#SAMtools: #"-S": Input format is SAM; "-h": Include header in output;
	#"-u": Output as uncompressed BAM; #"-F4": Do not output unmapped reads
	samtools view -ShuF4 $intermediateSAM | samtools sort - -o $sortedBAM

	#Create index file
	samtools index $sortedBAM
	
	#Select only reads that align to known regions of human and mouse genomes
	#"hg38"=human genome reference genome; "mm9"=mouse genome reference genome
	if [ $index == "hg38" ] || [ $index == "mm9" ]; then
		#De-duplicate reads by saving one per UMI
		umitools rmdup $sortedBAM temporary.bam > $BED
	
		#Create index file
		samtools index temporary.bam
		
		#Convert BAM to SAM file for processing
		samtools view -h -o temporary.sam temporary.bam
	
		#Remove reads that align to unidentified genomic regions
		#These unidentified regions include "EBV," "random," and "Un"
		sed '/chrEBV/d;/random/d;/chrUn/d' temporary.sam > filtered.sam
	
		#Convert SAM back to BAM file again
		samtools view -Sb filtered.sam > $finalBAM
	else
		#De-duplicate reads by saving one per UMI
		umitools rmdup $sortedBAM $finalBAM > $BED
	fi
	
	#Create index file
	samtools index $finalBAM
		
	#Remove intermediate and temporary files
	rm -f temporary.bam temporary.bam.bai temporary.sam \
	filtered.sam $sortedBAM $sortedBAM.bai $intermediateSAM
		
	#Notify user that alignment step is complete for which samples
	echo "Alignment of $sample to $index reference genome is complete"
done
