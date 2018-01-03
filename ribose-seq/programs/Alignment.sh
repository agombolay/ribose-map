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
		-s Sample name (e.g., FS100)
		-f Input Read 1 FASTQ filename
		-u UMI (e.g., NNNNNNNN or NNNNXXXNNNN)
		-b Barcode contained within UMI (e.g., TGA)
		-i Basename of Bowtie2 index (e.g., sacCer2)
		-m Minimum length of read to retain (e.g., 50)
		-d Local directory (e.g., /path/to/Ribose-Map)"
}

while getopts "u:m:i:f:s:b:d:h" opt; do
    	case "$opt" in
		u ) UMI=$OPTARG ;;
		m ) min=$OPTARG ;;
		i ) idx=$OPTARG ;;
		f ) read1=$OPTARG ;;
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
index=$directory/Indices/$idx
Fastq1=$directory/FASTQ-Files/$read1

#Output directory
output=$directory/Results/$idx/$sample/Alignment

#Create directory
mkdir -p $output

#############################################################################################################################
#Extract UMI from 5' ends of reads
umi_tools extract -I $Fastq1 -p $UMI -v 0 -S $output/UMI.fq

#Filter FASTQ file based on barcode sequence
grep --no-group-separator -B1 -A2 ^$barcode $output/UMI.fq > $output/filtered.fq

#Trim Illumina adapters and remove barcode from 5' end of reads
trim_galore --gzip --no_report_file --length $min --clip_R1 3 $output/filtered.fq -o $output

#############################################################################################################################
#Align reads to reference genome and save Bowtie2 statistics log file
bowtie2 -x $index -U $output/filtered_trimmed.fq.gz 2> $output/Align.log -S $output/mapped.sam
			
#Extract mapped reads, convert SAM file to BAM format, and sort and index BAM file
samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam && samtools index $output/sorted.bam
	
#############################################################################################################################		
#Remove PCR duplicates
umi_tools dedup -I $output/sorted.bam -v 0 | samtools sort - -o $output/$sample.bam && samtools index $output/$sample.bam
			
#############################################################################################################################
#Calculate percentage of reads that contain correct barcode sequence
x=$(echo $((`wc -l < $output/filtered.fq` / 4))/$((`wc -l < $output/UMI.fq` / 4)))

#Calculate percentage of reads that remain after de-duplication
y=$(echo "$(samtools view -c $output/$sample.bam)/$(samtools view -c $output/sorted.bam)")

#Save info about percentage of reads that contain correct barcode sequence
echo -e "Percentage of reads with barcode: $(echo "$x*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/Barcode.log

#Save info about percentage of reads that remain after de-duplication
echo -e "Percentage of reads that are unique: $(echo "$y*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/Unique.log
		
#############################################################################################################################
#Notify user alignment step is complete for input sample
echo "Trimming, alignment, and de-duplication of $sample is complete"

#Remove temporary files
rm -f $output/UMI.fq $output/filtered*.fq* $output/mapped.sam $output/sorted.bam*
