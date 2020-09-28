#!/usr/bin/env bash

#Author: Alli Lauren Gombolay
#Processes UMI and/or barcode if necessary, aligns reads to reference genome, and de-duplicates aligned reads

#############################################################################################################################
#Load config file
. "$1"

#Create output directory and remove any old files
output=$repository/results/$sample/alignment
rm -r $output; mkdir -p $output

#############################################################################################################################
if [[ ! $read2 ]]; then
	
	if [[ ! $pattern ]]; then

		if [[ ! $sort ]]; then
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -U $read1 -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam -o $output/$sample.bam

		elif [[ $sort ]]; then
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -U $read1 -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam	
		fi
		
	elif [[ $pattern ]]; then
		
		umi_tools extract -v 0 -I $read1 -p $pattern -S $output/extracted1.fq
		
		if [[ ! $barcode ]]; then
		
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -U $output/extracted1.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		
		elif [[ $barcode ]]; then
			
			grep -B 1 -A 2 ^$barcode $output/extracted1.fq | sed '/^--$/d' | seqtk trimfq -b ${#barcode} - > $output/demultiplexed1.fq
  
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -U $output/demultiplexed1.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
	
elif [[ $read2 ]]; then
	
	if [[ ! $pattern ]]; then

		if [[ ! $sort ]]; then
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -1 $read1 -2 $read2 -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam -o $output/$sample.bam

		elif [[ $sort ]]; then
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -1 $read1 -2 $read2 -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
		
	elif [[ $pattern ]]; then
		
		umi_tools extract -v 0 -I $read1 -p $pattern -S $output/extracted1.fq --read2-in=$read2 --read2-out=$output/extracted2.fq
  
		if [[ ! $barcode ]]; then
		
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -1 $output/extracted1.fq -2 $output/extracted2.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		
		elif [[ $barcode ]]; then
		
			grep -B 1 -A 2 ^$barcode $output/extracted1.fq | sed '/^--$/d' | seqtk trimfq -b ${#barcode} - > $output/demultiplexed1.fq
			
			bowtie2 -N $mismatches --seed $seed --threads $threads -x $basename -1 $output/demultiplexed1.fq -2 $output/extracted2.fq -S $output/aligned.sam 2> $output/alignment.log
			samtools view -bS -@ $threads $output/aligned.sam | samtools sort - -@ $threads -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -@ $threads -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
fi

#############################################################################################################################
#Print status
echo "Status: Alignment Module for $sample is complete"

if [[ $pattern ]]; then
	#Reads after de-duplication (%)
	duplication_statistics=$( bc -l <<< $(samtools view -c -F 4 $output/$sample.bam)/$(samtools view -c -F 4 $output/sorted.bam)*100 | xargs printf "%.*f\n" 2)
	echo "Reads after de-duplication (%)": $duplication_statistics >> $output/alignment.log
	
elif [[ $barcode ]]; then
	#Reads containing barcode (%)
	barcode_statistics=$( bc -l <<< $(wc -l $output/demultiplexed1.fq | awk '{print $1 / 4}')/$(wc -l $read1 | awk '{print $1 / 4}')*100 | xargs printf "%.*f\n" 2 )
	echo "Reads containing barcode (%)": $barcode_statistics >> $output/alignment.log
fi
