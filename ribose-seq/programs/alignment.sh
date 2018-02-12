#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Pre-processing (if UMI and/or barcode)
#2. Alignment or SE or PE reads to reference
#3. De-duplication based on UMI and chr coords

#Note: Input FASTQ files must be located in Ribose-Map 'fastqs' directory (Ribose-Map/fastqs)
#Note: Bowtie2 index files must be located in Ribose-Map 'indexes' directory (Ribose-Map/indexes)

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [options]
		Required:
		-d Ribose-Map repository
		-s Name of sequenced library
		-i Basename of Bowtie2 index
		-f Input Read 1 FASTQ filename
		-r Input Read 2 FASTQ filename
		Optional:
		-u UMI (e.g., NNNNNNNN or NNNNXXXNNNN)
		-b Molecular barcode in UMI (e.g., TGA)"
}

while getopts "u:i:f:r:s:b:d:h" opt; do
    	case "$opt" in
		u ) UMI=$OPTARG ;;
		i ) index=$OPTARG ;;
		f ) read1=$OPTARG ;;
		r ) read2=$OPTARG ;;
		s ) sample=$OPTARG ;;
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
index=$directory/indexes/$idx
fastq1=$directory/fastqs/$read1
fastq2=$directory/fastqs/$read2

#Create output directory
output=$directory/results/$sample/alignment

#Create directory and remove old files
mkdir -p $output; rm -f $output/*.{bam,bai,log}

#############################################################################################################################
if [[ ! $read2 ]]; then
	
	if [[ ! $umi ]]; then
		bowtie2 -x $idx -U $fastq1 -S $output/map.sam 2> $output/align.log
		
		samtools view -bS -F260 $output/map.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam	

	elif [[ $umi ]]; then
		umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/extract.fq
		
		if [[ ! $barcode ]]; then
			bowtie2 -x $idx -U $output/extract.fq -S $output/map.sam 2> $output/align.log

			samtools view -bS -F260 $output/map.sam | samtools sort - -o $output/sort.bam
			samtools index $output/sort.bam
	
			umi_tools dedup -v 0 -I $output/sort.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
		
		elif [[ $barcode ]]; then
			grep -B 1 -A 2 ^$barcode $output/extract.fq | sed '/^--$/d' \
			| awk 'NR%2 == 0 {sub(/^.{'${#barcode}'}/,"")} {print}' > $output/filter.fq
		
			bowtie2 -x $idx -U $output/filter.fq -S $output/map.sam 2> $output/align.log
			
			samtools view -bS -F260 $output/map.sam | samtools sort - -o $output/sort.bam
			samtools index $output/sort.bam
	
			umi_tools dedup -v 0 -I $output/sort.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
	
elif [[ $read2 ]]; then
	
	if [[ ! $umi ]]; then
		bowtie2 -x $idx -1 $fastq1 -2 $fastq2 -S $output/map.sam 2> $output/align.log
		
		samtools view -bS -f67 -F260 $output/map.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	
	elif [[ $umi ]]; then
	
		umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/extract1.fq --read2-in=$fastq2 \
		--read2-out=$output/extract2.fq
	
		if [[ ! $barcode ]]; then
			bowtie2 -x $idx -1 $output/extract1.fq -2 $output/extract2.fq -S $output/map.sam \
			2> $output/align.log
		
			samtools view -bS -f67 -F260 $output/map.sam | samtools sort - -o $output/sort.bam
			samtools index $output/sort.bam
	
			umi_tools dedup -v 0 --paired -I $output/sort.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
	
		elif [[ $barcode ]]; then
			grep -B 1 -A 2 ^$barcode $output/extract.fq | sed '/^--$/d' \
			| awk 'NR%2 == 0 {sub(/^.{'${#barcode}'}/,"")} {print}' > $output/filter.fq
		
			bowtie2 -x $idx -1 $output/filter.fq -2 $output/extract2.fq  -S $output/map.sam \
			2> $output/align.log
		
			samtools view -bS -f67 -F260 $output/map.sam | samtools sort - -o $output/sort.bam
			samtools index $output/sort.bam
	
			umi_tools dedup -v 0 --paired -I $output/sort.bam | samtools sort - -o $output/$sample.bam
			samtools index $output/$sample.bam
		fi
	fi
fi
		
#############################################################################################################################
#Print completion status
echo "Status: Program complete for $sample"

#Remove temporary files from directory
rm -f $output/*.{fq,fq.gz,sam} $output/sorted.{bam,bai}
