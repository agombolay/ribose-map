#!/usr/bin/env bash

#Author: Alli Gombolay
#This program raligns trimmed reads to reference genome using Bowtie2, and de-duplicates reads based on UMI's

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 1_Alignment.sh [-i] 'FASTQ' [-b] 'Index' [-d] 'Directory' [-h]
		-i Sample names (FS1.trimmed.v1 FS2.trimmed.v1 FS3.trimmed.v1 etc.) 
		-b Basename of Bowtie2 index to be searched (sacCer2, chrM, ecoli, hg38, etc.)
		-d Local directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Use getopts function to create the command-line options ([-i], [-b], [-d], and [-h])
while getopts "i:b:d:h" opt; do
    case "$opt" in
        #Specify input as arrays to allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	b ) index=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#############################################################################################################################
#Align reads to reference genome
for sample in ${sample[@]}; do
	
	#Input FASTQ files
	fastq=$directory/Sequencing-Results/$sample.fastq
	
	#Output directory
	output=$directory/ribose-seq/results/$index/$sample/Alignment/

	#Create directory if not present
    	mkdir -p $output

	#Intermediate files
	umiTrimmed=$output/$sample.UMI-trimmed.fastq.gz; intermediateSAM=$output/$sample.intermediate.sam;
	sortedBAM=$output/$sample.sorted.bam; reverseComplement=$output/$sample.reverse-complement.fastq

	#UMI length
	UMI=NNNNNNNN
	
	#BED file
	BED=$output/$sample.bed.gz
	
	#Final BAM file
	finalBAM=$output/$sample.bam

	#File containing Bowtie2 alignment statistics
	statistics=$output/$sample.Alignment-Statistics.txt

#############################################################################################################################
	#1. Reverse complement reads
	seqtk seq -r $fastq > $reverseComplement

	#2. Trim UMI from 3' ends of reads; compress file
	umitools trim --end 3 $reverseComplement $UMI | gzip -c > $umiTrimmed
	
	#Alternative Version: Trim UMI from 5' ends of reads
	#umitools trim --end 5 $fastq $UMI | gzip -c > $umiTrimmed

	#3. Align reads to reference genome using Bowtie2
	#Bash: "-": standard input; "2>": Redirect standard error; "-x": Index;
	#"-U": Unpaired input reads; "-S": Print alignment results in SAM format
	zcat $umiTrimmed | bowtie2 -x $index -U - -S $intermediateSAM 2> $statistics

	#Alternative Version: Align reads using Bowtie version 1
	#Bash: "-": standard input; "2>": Redirect standard error; "1>": Redirect standard output
	#Bowtie: "-m 1": Return only unique reads; "--sam": Print alignment results in SAM format
	#zcat $umiTrimmed | bowtie -m 1 --sam $index - 2> $statistics 1> $intermediateSAM
	
	#4. Convert SAM file to BAM and sort intermediate BAM file
	#SAMtools: #"-S": Input format is SAM; "-h": Include header in output;
	#"-u": Output as uncompressed BAM; #"-F4": Do not output unmapped reads
	samtools view -ShuF4 $intermediateSAM | samtools sort - $sortedBAM

	#5. Create index file
	samtools index $sortedBAM
	
	#Select only reads that align to known regions of human and mouse genomes
	#"hg38"=human genome reference genome; "mm9"=mouse genome reference genome
	if [ $index == "hg38" ] || [ $index == "mm9" ]; then
		#6. De-duplicate reads based on UMIs; compress file
		umitools rmdup $sortedBAM temporary.bam | gzip -c > $BED
	
		#7. Create index file
		samtools index temporary.bam
		
		#8. Convert BAM to SAM file for processing
		samtools view -h -o temporary.sam temporary.bam
	
		#9. Remove reads that align to unidentified genomic regions
		#These unidentified regions include "EBV," "random," and "Un"
		sed '/chrEBV/d;/random/d;/chrUn/d' temporary.sam > filtered.sam
	
		#10. Convert SAM back to BAM file again
		samtools view -Sb filtered.sam > $finalBAM
	
		#11. Create index file
		samtools index $finalBAM
	else
		#6. De-duplicate reads based on UMIs; compress file
		umitools rmdup $sortedBAM $finalBAM | gzip -c > $BED
	
		#7. Create index file
		samtools index $finalBAM
	fi
	
	#Remove intermediate and temporary files
	rm -f temporary.bam temporary.bam.bai temporary.sam \
	filtered.sam $sortedBAM $sortedBAM.bai $intermediateSAM
		
	#Notify user that alignment step is complete for which samples
	echo "Alignment of $sample to $index reference genome is complete"
done
