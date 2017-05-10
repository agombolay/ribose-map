#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program raligns trimmed reads to reference genome using Bowtie2 and de-duplicates reads based on UMI's
#Note: FASTQ files must be located in the users's Sequencing-Results folder (/LocalDirectory/Sequencing-Results)

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [-i] 'Sample(s)' [-u] 'UMI' [-m] 'Min' [-p] 'Path' [-b] 'Index' [-d] 'Directory' [-h]
		-i Sample name(s) (FS1, FS2, FS3 etc.)
		-u Length of UMI (NNNNNNNN or NNNNNNNNNNN)
		-m Minimum length of read to retain after trimming
		-p /projects/home/agombolay3/data/bin/Trimmomatic-0.36
		-b Basename of Bowtie2 index (sacCer2, pombe, ecoli, mm9, hg38, etc.)
		-d Directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
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
	
	#Output directory
	output=$directory/ribose-seq/results/$index/$sample/Alignment

	#Create directory
    	mkdir -p $output

	#Intermediate files
	umiTrimmed=$output/$sample.UMI-trimmed.fastq.gz; intermediateSAM=$output/$sample.intermediate.sam;
	sortedBAM=$output/$sample.sorted.bam; reverseComplement=$output/$sample.reverse-complement.fastq;
	mappedBAM=$output/$sample.mapped.bam; unmappedBAM=$output/$sample.unmapped.bam;
	unmappedFASTQ=$output/$sample.unmapped.fastq;
	
	#BED file
	BED=$output/$sample.bed.gz
	
	#Final BAM file
	finalBAM=$output/$sample.bam

	#File containing Bowtie2 alignment statistics
	statistics=$output/$sample.Alignment-Statistics.txt

#############################################################################################################################
	#Trim FASTQ files based on quality and Illumina adapter content
	java -jar $path/trimmomatic-0.36.jar SE -phred33 $fastq $output/$sample-trimmed.fastq \
	ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:10 SLIDINGWINDOW:5:15 MINLEN:$MIN

	#Trim UMI from 5' ends of reads (add UMI into read name for further processing)
	umitools trim --end 5 $output/$sample-trimmed.fastq $UMI | gzip -c > $umiTrimmed

	#Reverse complement reads
	zcat $umiTrimmed | seqtk seq -r - > $reverseComplement
	
	#Align reads to reference genome using Bowtie
	#zcat $umiTrimmed | bowtie -m 1 $index - -S $intermediateSAM 2> $statistics
	bowtie2 -x $index -U $reverseComplement -S $intermediateSAM 2> $statistics
	
	#Convert SAM file to BAM and sort intermediate BAM file
	#SAMtools: #"-S": Input format is SAM; "-h": Include header in output;
	#"-u": Output as uncompressed BAM; #"-F4": Do not output unmapped reads
	#samtools view -ShuF4 $intermediateSAM | samtools sort - -o $sortedBAM
	samtools view -Shu $intermediateSAM | samtools sort - -o $sortedBAM

	#Save mapped reads to BAM file
	samtools view -bF4 $sortedBAM > $mappedBAM
	
	#Save unmapped reads to FASTQ file
	samtools view -bf4 $sortedBAM > $unmappedBAM
	bamToFastq -i $unmappedBAM -fq $unmappedFASTQ
	
	#Create index file
	samtools index $mappedBAM
	
	#De-duplicate reads by saving one per UMI
	umitools rmdup $mappedBAM $finalBAM > $BED
	
	#Create index file
	samtools index $finalBAM
		
	#Remove intermediate and temporary files
	rm -f $sortedBAM $sortedBAM.bai $intermediateSAM
		
	#Notify user that alignment step is complete for which samples
	echo "Alignment of $sample to $index reference genome is complete"
done
