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
		-d Ribose-Map directory
		-s Name of sequenced library
		-i Basename of Bowtie2 index
		-f Input Read 1 FASTQ filename
		-r Input Read 2 FASTQ filename
		Optional:
		-u UMI (e.g., NNNNNNNN or NNNNXXXNNNN)
		-b Molecular barcode in UMI (e.g., TGA)"
}

while getopts "h:u:i:f:r:s:b:d" opt; do
    	case "$opt" in
		h ) usage ;;
		u ) UMI=$OPTARG ;;
		i ) index=$OPTARG ;;
		f ) read1=$OPTARG ;;
		r ) read2=$OPTARG ;;
		s ) sample=$OPTARG ;;
		b ) barcode=$OPTARG ;;
		d ) directory=$OPTARG ;;
    	esac
done

#############################################################################################################################
#Input files
fastq1=$directory/fastqs/$read1
fastq2=$directory/fastqs/$read2
index=$directory/references/$idx

#Create output directory
output=$directory/results/$sample/alignment

#Create directory and remove old files
mkdir -p $output; rm -f $output/*.{bam,bai,log}

#############################################################################################################################
if [[ ! $read2 ]]; then
	
	if [[ ! $umi ]]; then
		bowtie2 -x $idx -U $fastq1 -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam	

	elif [[ $umi ]]; then
		bowtie2 -x $idx -U $fastq1 -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
		samtools index $output/sort.bam
	
		umi_tools dedup -v 0 -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	fi
	
elif [[ $read2 ]]; then
	
	if [[ ! $umi ]]; then
		bowtie2 -x $idx -1 $fastq1 -2 $fastq2 -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -f67 -F260 $output/mapped.sam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	
	elif [[ $umi ]]; then	
		bowtie2 -x $idx -1 $fastq1 -2 $fastq2 -S $output/mapped.sam 2> $output/alignment.log
		samtools view -bS -f67 -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
		samtools index $output/sorted.bam
	
		umi_tools dedup -v 0 --paired -I $output/sorted.bam | samtools sort - -o $output/$sample.bam
		samtools index $output/$sample.bam
	fi
fi

#Calculate % of reads that remain after de-duplication step
y=$(echo "$(samtools view -c $output/$sample.bam)/$(samtools view -c $output/sorted.bam)")

#Save info about % of reads that remain after de-duplication step
echo -e "Percentage of reads that are unique: $(echo "$y*100" | bc -l)%" > $output/unique.log

#############################################################################################################################
#Print completion status
echo "Status: Program complete for $sample"

#Remove temporary files from directory
rm -f $output/*.{fq,sam} $output/sorted.{bam,bai}
