#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Pre-processing (if UMI and/or barcode)
#2. Alignment or SE or PE reads to reference
#3. De-duplication based on UMI and chr coords

#############################################################################################################################
#Load config file
. "$1"

#Create output directory
output=$repository/results/$sample/alignment; rm -rf $output; mkdir -p $output

#############################################################################################################################
if [[ ! $read2 ]]; then
	
	if [[ ! $pattern ]]; then
		
		if [[ ! $sort ]]; then
			bowtie2 --threads $threads -x $basename -U $read1 -S $output/aligned.sam 2> $output/alignment.log
			samtools view -b -S -@ $threads $output/aligned.sam -o $output/$sample.bam
		
		elif [[ $sort ]]; then
			bowtie2 --threads $threads -x $basename -U $read1 -S $output/aligned.sam 2> $output/alignment.log
			samtools view -b -S -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam	
		
	elif [[ $pattern ]]; then
		
		umi_tools extract -v 0 -I $read1 -p $pattern -S $output/extracted1.fq
		
		if [[ ! $barcode ]]; then
		
			bowtie2 --threads $threads -x $basename -U $output/extracted1.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -b -S -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		
		elif [[ $barcode ]]; then
			
			grep -B 1 -A 2 ^$barcode $output/extracted1.fq | sed '/^--$/d' | seqtk trimfq -b ${#barcode} - > $output/demultiplexed1.fq
  
			bowtie2 --threads $threads -x $basename -U $output/demultiplexed1.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -b -S -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
	
elif [[ $read2 ]]; then
	
	if [[ ! $pattern ]]; then
	
		bowtie2 --threads $threads -x $basename -1 $read1 -2 $read2 -S $output/aligned.sam 2> $output/alignment.log
		samtools view -b -S -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/$sample.bam
		samtools index $output/$sample.bam
	
	elif [[ $pattern ]]; then
		
		umi_tools extract -v 0 -I $read1 -p $pattern -S $output/extracted1.fq --read2-in=$read2 --read2-out=$output/extracted2.fq
  
		if [[ ! $barcode ]]; then
		
			bowtie2 --threads $threads -x $basename -1 $output/extracted1.fq -2 $output/extracted2.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -b -S -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		
		elif [[ $barcode ]]; then
		
			grep -B 1 -A 2 ^$barcode $output/extracted1.fq | sed '/^--$/d' | seqtk trimfq -b ${#barcode} - > $output/demultiplexed1.fq
			
			bowtie2 --threads $threads -x $basename -1 $output/demultiplexed1.fq -2 $output/extracted2.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -b -S -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
fi

#############################################################################################################################
#Print status
echo "Status: Alignment Module for $sample is complete"
