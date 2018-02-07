#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Pre-processing, 2. Alignment, and 3. De-duplication based on UMI

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
		-a Custom adapter sequence to remove
		-m Minimum length of reads to retain
		-u UMI (e.g., NNNNNNNN or NNNNXXXNNNN)
		-b Molecular barcode in UMI (e.g., TGA)"
}

while getopts "u:m:i:f:s:a:b:d:h" opt; do
    	case "$opt" in
		u ) UMI=$OPTARG ;;
		m ) min=$OPTARG ;;
		i ) idx=$OPTARG ;;
		f ) read1=$OPTARG ;;
		s ) sample=$OPTARG ;;
		a ) adapter=$OPTARG ;;
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
index=$directory/indexes/$idx;
fastq1=$directory/fastqs/$read1; fastq2=$directory/fastqs/$read2

#Create output directory and remove old directory if present
output=$directory/results/$sample/alignment; mkdir -p $output; rm -rf $output/*

#############################################################################################################################
#Extract UMI from 5' ends of reads
if [[ $umi ]] && [[ ! $read2 ]]; then
	umi_tools extract -v 0 -I $fastq1 --bc-pattern=$UMI -S processed.fq.gz
	
elif [[ $umi ]] && [[ $read2 ]]; then
	umi_tools extract -v 0 -I $fastq1 --bc-pattern=$UMI --read2-in=$fastq2 \
	--stdout=processed1.fq.gz --read2-out=processed2.fq.gz
fi

#############################################################################################################################
if [[ $barcode ]]; then
	#Filter FASTQ file based on barcode sequence
	grep --no-group-separator -B1 -A2 ^$barcode $output/UMI.fq > $output/filtered.fq
fi

#############################################################################################################################
if [[ ! $adapter ]] && [[ $illumina ]]; then
	#Trim Illumina and remove barcode from 5' end of reads
	trim_galore --gzip --length $min --clip_R1 3 $output/filtered.fq -o $output

elif [[ $adapter ]] && [[ $illumina ]]; then
	#Trim Illumina/custom adapters and remove barcode from 5' end of reads
	trim_galore --gzip --length $min --clip_R1 3 -a $adapter $output/filtered.fq -o $output
	
	trim_galore --gzip --length $min --clip_R1 3 $output/filtered.fq -o $output
	trim_galore --gzip --paired --length $min --clip_R1 3 -a $adapter $output/filtered.fq -o $output
fi

#############################################################################################################################
if [[ $type == 'se' ]]; then
	#Align reads to reference genome and save Bowtie2 statistics log file
	bowtie2 -x $index -U $output/filtered_trimmed.fq.gz 2> $output/alignment.log -S $output/mapped.sam
	
	#Extract mapped reads, convert SAM file to BAM format, and sort/index BAM file
	samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam && samtools index $output/sorted.bam

elif [[ $type == 'pe' ]]; then
	#Align reads to reference genome and save Bowtie2 statistics log file
	bowtie2 -x $index -1 $output/.fq.gz -2 $output/.fq.gz 2> $output/alignment.log -S $output/mapped.sam

fi

#############################################################################################################################
if [[ $umi ]] && [[ ! $read2 ]]; then
	#Remove PCR duplicates based on UMI and mapping coordinates and sort/index BAM file
	umi_tools dedup -I $output/sorted.bam -v 0 | samtools sort - -o $output/$sample.bam && samtools index $output/$sample.bam

elif [[ $umi ]] && [[ $read2 ]]; then
	umi_tools dedup --paired -I mapped.bam -S deduplicated.bam
fi

#############################################################################################################################
#Calculate percentage of reads that contain correct barcode sequence
x=$(echo $((`wc -l < $output/filtered.fq` / 4))/$((`wc -l < $output/UMI.fq` / 4)))

#Calculate percentage of reads that remain after de-duplication
y=$(echo "$(samtools view -c $output/$sample.bam)/$(samtools view -c $output/sorted.bam)")

#Save info about percentage of reads that contain correct barcode sequence
echo -e "Percentage of reads with barcode: $(echo "$x*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/barcode.log

#Save info about percentage of reads that remain after de-duplication
echo -e "Percentage of reads that are unique: $(echo "$y*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/unique.log
		
#############################################################################################################################
#Print completion status
echo "Status: Program complete for $sample"

#Remove temporary files from directory
rm -f $output/*.{fq,fq.gz,sam} $output/sorted.{bam,bai}
