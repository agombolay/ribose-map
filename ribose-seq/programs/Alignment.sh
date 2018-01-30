#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Pre-process reads, 2. align reads to reference genome, and 3. de-duplicate reads based on UMI

#Note: Input FASTQ files must be located in Ribose-Map 'fastqs' directory (filepath/Ribose-Map/fastqs)
#Note: Bowtie2 index files must be located in Ribose-Map 'indexes' directory (filepath/Ribose-Map/indexes)

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [options]
		-s Sample name (e.g., FS100)
		-f Input Read 1 FASTQ filename
		-r Input Read 2 FASTQ filename
		-u UMI (e.g., NNNNNNNN or NNNNXXXNNNN)
		-b Barcode contained within UMI (e.g., TGA)
		-i Basename of Bowtie2 index (e.g., sacCer2)
		-a Extra adapter sequence to be removed from reads
		-m Minimum length of reads to retain after trimming
		-d Ribose-Map directory (e.g., /path/to/Ribose-Map)"
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
index=$directory/indexes/$idx
fastq1=$directory/fastqs/$read1
fastq2=$directory/fastqs/$read2

#Output directory
output=$directory/results/$idx/$sample/alignment

#############################################################################################################################
#Create directory
mkdir -p $output

#Remove any old files
rm -f $output/*.{bam,log}

#############################################################################################################################
#Extract UMI from 5' ends of reads
umi_tools extract -I $fastq1 -p $UMI -v 0 -S $output/UMI.fq

#Filter FASTQ file based on barcode sequence
grep --no-group-separator -B1 -A2 ^$barcode $output/UMI.fq > $output/filtered.fq

if [[ ! $adapter ]]; then
	#Trim Illumina and remove barcode from 5' end of reads
	trim_galore --gzip --length $min --clip_R1 3 $output/filtered.fq -o $output

elif [[ $adapter ]]; then
	#Trim Illumina/custom adapters and remove barcode from 5' end of reads
	trim_galore --gzip --length $min --clip_R1 3 -a $adapter $output/filtered.fq -o $output
fi

#############################################################################################################################
#Align reads to reference genome and save Bowtie2 statistics log file
bowtie2 -x $index -U $output/filtered_trimmed.fq.gz 2> $output/alignment.log -S $output/mapped.sam
			
#Extract mapped reads, convert SAM file to BAM format, and sort/index BAM file
samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam && samtools index $output/sorted.bam
	
#Remove PCR duplicates based on UMI and mapping coordinates and sort/index BAM file
umi_tools dedup -I $output/sorted.bam -v 0 | samtools sort - -o $output/$sample.bam && samtools index $output/$sample.bam
			
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
#Notify user alignment step is complete for input sample
echo "Trimming, alignment, and de-duplication of $sample is complete"

#Remove temporary files
rm -f $output/*.fq $output/*.fq.gz $output/mapped.sam $output/sorted.bam*
