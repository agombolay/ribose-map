#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program aligns trimmed reads to reference genome using Bowtie2 and de-duplicates reads based on UMI
#Note1: FASTQ files must be located in users's FASTQ-Files folder (/LocalDirectory/Ribose-Map/FASTQ-Files)
#Note2: rNMP is the reverse complement of the 5' base of the sequenced read in FASTQ file

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [options]
		-s Sample name(s) (e.g., FS1, FS2, FS3)
		-f Input Read 1 FASTQ filename (forward)
		-r Input Read 2 FASTQ filename (reverse)
		-u Length of UMI (e.g., NNNNNNNN or NNNNNNNNNNN)
		-b Barcode contained within UMI (e.g., ..TGA......)
		-m Minimum length of read to retain after trimming (e.g., 50)
		-p Path (e.g., /projects/home/agombolay3/data/bin/Trimmomatic-0.36)
		-t Type of Illumina Sequencing (e.g., SE = Single end, PE = Paired end)
		-i Basename of Bowtie2 index (e.g., sacCer2, pombe, ecoli, mm9, hg38, etc.)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository/Ribose-Map)"
}

while getopts "s:u:p:t:i:f:r:m:b:d:h" opt; do
    	case "$opt" in
        	#Allow multiple input arguments
        	s ) sample=($OPTARG) ;;
		#Allow only one input argument
		u ) UMI=$OPTARG ;;
		p ) path=$OPTARG ;;
		t ) type=$OPTARG ;;
		i ) index=$OPTARG ;;
		f ) read1=$OPTARG ;;
		r ) read2=$OPTARG ;;
		m ) minimum=$OPTARG ;;
		b ) barcode=$OPTARG ;;
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
#Input files
Fastq1=$directory/FASTQ-Files/$read1
Fastq2=$directory/FASTQ-Files/$read2

#Output directory
output=$directory/Results/$index/$sample/Alignment

#Create directory
mkdir -p $output

#############################################################################################################################
for sample in ${sample[@]}; do

	#Remove old files
	rm -f $output/$sample.bam*

	#Single End Reads
	if [[ $type == "SE" ]]; then
	
		#Reverse complement reads
		cat $Fastq1 | seqtk seq -r - > $output/Reverse.fq
	
		#Extract UMI from 3' ends of reads and append to read name
		umi_tools extract -I $output/Reverse.fq -p $UMI --3prime -v 0 -S $output/Extract.fq
	
		#Trim/drop reads based on quality, adapter content, and length
		java -jar $path/trimmomatic-0.36.jar SE -trimlog $output/Trimmomatic.log $output/Extract.fq \
		$output/Read1.fq ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 LEADING:10 MINLEN:$minimum
		
		#Align reads to reference genome and save Bowtie2 statistics log file
		bowtie2 -x $index -U $output/Read1.fq 2> $output/Bowtie2.log > $output/mapped.sam
			
		#Extract mapped reads, convert SAM file to BAM format, and sort BAM file
		samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
			
		#Index BAM file
		samtools index $output/sorted.bam
			
		if [[ -n $UMI ]] && [[ -z $barcode ]]; then
			
			#Remove PCR duplicates based on UMI and genomic start position and sort BAM file
			umi_tools dedup -I $output/sorted.bam -v 0 | samtools sort - -o $output/$sample.bam
			
			#Index BAM file
			samtools index $output/$sample.bam
		
		elif [[ -n $UMI ]] && [[ -n $barcode ]]; then
			
			#Remove PCR duplicates based on UMI and genomic start position and sort BAM file
			umi_tools dedup -I $output/sorted.bam -v 0 | samtools sort - -o $output/deduped.bam
			
			#Filter BAM file based on barcode
			samtools view -h $output/deduped.bam -o $output/deduped.sam
			grep -e "_$barcode" -e '^@' $output/deduped.sam > $output/filtered.sam
			samtools view $output/filtered.sam -bS | samtools sort -o $output/$sample.bam
			
			#Index BAM file
			samtools index $output/$sample.bam
		
		fi
	fi

	#Paired End Reads
	if [[ $type == "PE" ]]; then
	
		#Reverse complement reads
		cat $Fastq1.fq | seqtk seq -r - > $output/Read1.fq
		cat $Fastq2.fq | seqtk seq -r - > $output/Read2.fq
	
		#Extract UMI from 3' ends of reads and append to read name
		umi_tools extract -I $output/Read1.fq -p $UMI --3prime -v 0 -S $output/Extract.fq
	
		#Trim/drop reads based on quality, adapter content, and length
		java -jar $path/trimmomatic-0.36.jar PE -trimlog $output/Trimmomatic.log $output/Extract.fq \
		$output/Read2.fq $output/Paired1.fq $output/Unpaired1.fq $output/Paired2.fq $output/Unpaired2.fq \
		ILLUMINACLIP:$path/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 MINLEN:$minimum
		
		#Align reads to reference genome and save Bowtie2 statistics log file
		bowtie2 -x $index -1 $output/Paired1.fq -2 $output/Paired2.fq -S $output/mapped.sam
		
		#Extract mapped reads, convert SAM file to BAM format, and sort BAM file
		samtools view -bS -f66 -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
		
		#Index BAM file
		samtools index $output/sorted.bam
		
		if [[ -n $UMI ]] && [[ -z $barcode ]]; then
		
			#Remove PCR duplicates based on UMI and genomic start position and sort BAM file
			umi_tools dedup -I $output/sorted.bam -v 0 | samtools sort - -o $output/$sample.bam
		
			#Index BAM file
			samtools index $output/$sample.bam
	
		elif [[ -n $UMI ]] && [[ -n $barcode ]]; then
		
			#Remove PCR duplicates based on UMI and genomic start position and sort BAM file
			umi_tools dedup -I $output/sorted.bam -v 0 | samtools sort - -o $output/deduped.bam
		
			#Filter BAM file based on barcode
			samtools view -h $output/deduped.bam -o $output/deduped.sam
			grep -e "_$barcode" -e '^@' $output/deduped.sam > $output/filtered.sam
			samtools view $output/filtered.sam -bS | samtools sort -o $output/$sample.bam
		
			#Index BAM file
			samtools index $output/$sample.bam
	
		fi
	fi

#############################################################################################################################
	#Notify user alignment step is complete for input sample
	echo "Trimming, alignment, and de-duplication of $sample is complete"

	#Remove temporary files
	rm -f $output/Reverse.fq $output/Extract.fq $output/Read1.fq $output/Read2.fq \
	$output/mapped.sam $output/sorted.bam* $output/deduped.bam* $output/filtered.sam \
	$output/Paired1.fq $output/Unpaired1.fq $output/Paired2.fq $output/Unpaired2.fq

done
