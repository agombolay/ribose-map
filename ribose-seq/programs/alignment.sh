#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Pre-processing (if UMI and/or barcode)
#2. Alignment or SE or PE reads to reference
#3. De-duplication based on UMI and chr coords

#Note: Input FASTQ files must be located in Ribose-Map 'fastqs' directory (Ribose-Map/fastqs)
#Note: Bowtie2 index files must be located in Ribose-Map 'indexes' directory (Ribose-Map/indexes)

. /data2/users/agombolay3/Ribose-Map/config.txt

output=$directory/results/$sample/alignment
mkdir -p $output

#############################################################################################################################
if [[ ! $read2 ]]; then
	
	if [[ ! $umi ]]; then
		bowtie2 -x $index -U $read1_fastq -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam	

	elif [[ $umi ]]; then
		bowtie2 -x $index -U $read1_fastq -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
		samtools index $output/sort.bam
	
		umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	fi
	
elif [[ $read2 ]]; then
	
	if [[ ! $umi ]]; then
		bowtie2 -x $index -1 $read1_fastq -2 $read2_fastq -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -f67 -F260 $output/mapped.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	
	elif [[ $umi ]]; then	
		bowtie2 -x $index -1 $read1_fastq -2 $read2_fastq -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -f67 -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
		samtools index $output/sorted.bam
	
		umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	fi
fi

#############################################################################################################################
#Print completion status
#echo "Status: Program complete for $sample"

#Remove temporary files from directory
#rm -f $output/*.{fq,sam} $output/sorted.{bam,bai}
