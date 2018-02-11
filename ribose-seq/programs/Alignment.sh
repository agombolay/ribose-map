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
index=$directory/indexes/$idx;
fastq1=$directory/fastqs/$read1
fastq2=$directory/fastqs/$read2

#Create output directory
output=$directory/results/$sample/alignment

#Create directory and remove old files
mkdir -p $output; rm -f $output/*.{bam,bai,log}

#############################################################################################################################
#Extract UMI from reads and append UMI to read name
if [[ $umi ]] && [[ ! $read2 ]]; then
	umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/processed1.fq.gz

elif [[ $umi ]] && [[ $read2 ]]; then
	umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/processed1.fq.gz \
	--read2-in=$fastq2 --read2-out=$output/processed2.fq.gz
fi

#Filter reads by barcode and remove barcode from reads
if [[ $barcode ]]; then
	grep --no-group-separator -B1 -A2 ^$barcode $output/processed1.fq.gz \ 
	| cutadapt - --cut ${#barcode} -o $output/output.fq.gz
fi

#############################################################################################################################
#Align reads to reference genome and save Bowtie2 log file
#Extract mapped reads, convert SAM file to BAM, and sort/index
if [[ ! $read2 ]]; then
	bowtie2 -x $index -U $output/filtered_trimmed.fq.gz 2> $output/alignment.log -S $output/mapped.sam
	samtools view -bS -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
	samtools index $output/sorted.bam
	
elif [[ $read2 ]]; then
	bowtie2 -x $index -1 $output/.fq.gz -2 $output/.fq.gz 2> $output/alignment.log -S $output/mapped.sam
	samtools view -bS -f67 -F260 $output/mapped.sam | samtools sort - -o $output/sorted.bam
	samtools index $output/sorted.bam
fi

#############################################################################################################################
#De-duplicate aligned reads based on UMI and chromosome coordinates
if [[ $umi ]] && [[ ! $read2 ]]; then
	umi_tools dedup -v 0 -I $output/mapped.bam | samtools sort - -o $output/mapped.bam

elif [[ $umi ]] && [[ $read2 ]]; then
	umi_tools dedup -v 0 --paired -I $output/mapped.bam | samtools sort - -o $output/mapped.bam
fi

#############################################################################################################################
#Calculate % of reads that contain correct barcode sequence
x=$(echo $((`wc -l < $output/filtered.fq` / 4))/$((`wc -l < $output/UMI.fq` / 4)))

#Calculate % of reads that remain after de-duplication step
y=$(echo "$(samtools view -c $output/$sample.bam)/$(samtools view -c $output/sorted.bam)")

#Save info about %of reads that contain correct barcode sequence
echo -e "Percentage of reads with barcode: $(echo "$x*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/barcode.log

#Save info about % of reads that remain after de-duplication step
echo -e "Percentage of reads that are unique: $(echo "$y*100" | bc -l | xargs printf "%.*f\n" 2)%" > $output/unique.log
		
#############################################################################################################################
#Print completion status
echo "Status: Program complete for $sample"

#Remove temporary files from directory
rm -f $output/*.{fq,fq.gz,sam} $output/sorted.{bam,bai}
