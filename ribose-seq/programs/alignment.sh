#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Pre-processing (if UMI and/or barcode)
#2. Alignment or SE or PE reads to reference
#3. De-duplication based on UMI and chr coords

. /data2/users/agombolay3/Ribose-Map/config.txt

output=$directory/results/$sample/alignment
mkdir -p $output

#############################################################################################################################
if [[ ! $read2 ]]; then
	
	if [[ ! $umi ]]; then
	
		bowtie2 -x $idx -U $read1 -S $output/aligned.sam 2> $output/alignment.log
		samtools view -bS -F260 $output/aligned.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam	

	elif [[ $umi ]]; then
		
		umi_tools extract -v 0 -I $read1 -p $UMI -S $output/umi.fq
		
		if [[ ! $barcode ]]; then
		
			bowtie2 -x $idx -U $output/umi.fq -S $output/aligned.sam 2> $output/align.log
			samtools view -bS -F260 $output/aligned.sam | samtools sort - -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
		
		elif [[ $barcode ]]; then
		
			grep -B 1 -A 2 --ignore-case ^$barcode $output/umi.fq | sed '/^--$/d' > $output/trim.fq
			awk 'NR % 2 == 0 {sub(/^.{'${#barcode}'}/,"")} {print}' $output/trim.fq > $output/filter.fq
  
			bowtie2 -x $idx -U $output/filter.fq -S $output/aligned.sam 2> $output/align.log
			samtools view -bS -F260 $output/aligned.sam | samtools sort - -o $output/sorted.bam
			samtools index $output/sort.bam
	
			umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
	
elif [[ $read2 ]]; then
	
	if [[ ! $umi ]]; then
	
		bowtie2 -x $idx -1 $read1 -2 $read2 -S $output/aligned.sam 2> $output/align.log
		samtools view -bS -f67 -F260 $output/aligned.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	
	elif [[ $umi ]]; then
		
		umi_tools extract -v 0 -I $read1 -p $UMI -S $output/umi1.fq --read2-in=$read2 --read2-out=$output/umi2.fq
  
		if [[ ! $barcode ]]; then
		
			bowtie2 -x $idx -1 $output/filter.fq -2 $output/umi2.fq -S $output/aligned.sam 2> $output/align.log
			
			samtools view -bS -f67 -F260 $output/aligned.sam | samtools sort - -o $output/sorted.bam
			samtools index $output/sorted.bam
	
			umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
		
		elif [[ $barcode ]]; then
		
			grep -B 1 -A 2 --ignore-case ^$barcode $output/umi1.fq | sed '/^--$/d' > $output/trim.fq
			awk 'NR % 2 == 0 {sub(/^.{'${#barcode}'}/,"")} {print}' $output/trim.fq > $output/filter.fq
			
			bowtie2 -x $idx -1 $output/filter.fq -2 $output/umi2.fq -S $output/aligned.sam 2> $output/align.log
			
			samtools view -bS -f67 -F260 $output/aligned.sam | samtools sort - -o $output/sorted.bam
			samtools index $output/sort.bam
	
			umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
fi

#############################################################################################################################
#Calculate % of reads that contain correct barcode sequence
#x=$(echo $(bc -l <<< "$(wc -l < $output/barcode.fq)/4")/$(bc -l <<< "$(wc -l < $output/umi_extracted1.fq)/4"))

#Save info about % of reads that remain after de-duplication step
#echo -e "Reads that are unique molecules: $(echo "$y*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/duplication.log

#Calculate % of reads that remain after de-duplication step
#y=$(echo $(bc -l <<< "$(samtools view -c < $output/$sample.bam)")/$(bc -l <<< "$(samtools view -c < $output/sorted.bam)"))

#Save info about % of reads that contain correct barcode sequence
#echo -e "Reads with molecular barcode, $barcode: $(echo "$x*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/barcode.log

#Print completion status
#echo "Status: Program complete for $sample"

#Remove temporary files from directory
#rm -f $output/*.{fq,sam} $output/sorted.{bam,bai}
