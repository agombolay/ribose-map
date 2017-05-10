#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program raligns trimmed reads to reference genome using Bowtie2 and de-duplicates reads based on UMI's
#Note: FASTQ files must be located in the users's Sequencing-Results folder (/LocalDirectory/Sequencing-Results)

#Usage statement
function usage () {
	echo "Usage: Alignment.sh [-i] 'Sample(s)' [-u] 'UMI' [-m] 'Min' [-p] 'Path' [-b] 'Index' [-d] 'Directory' [-h]
		-i Input sample(s) (e.g., FS1, FS2, FS3)
		-u Length of UMI (e.g., NNNNNNNN or NNNNNNNNNNN)
		-m Minimum length of read to retain after trimming (e.g., 50)
		-p Path (e.g., /projects/home/agombolay3/data/bin/Trimmomatic-0.36)
		-b Basename of Bowtie2 index (e.g., sacCer2, pombe, ecoli, mm9, or hg38)
		-d Directory (e.g., /projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Command-line options
while getopts "i:u:m:p:b:d:h" opt; do
    case "$opt" in
        #Allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Allow only one input argument
	u ) UMI=$OPTARG ;;
	m ) MIN=$OPTARG ;;
	p ) path=$OPTARG ;;
	b ) index=$OPTARG ;;
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
#Align reads to reference
for sample in ${sample[@]}; do
	
	#Input FASTQ files
	fastq=$directory/Sequencing-Results/$sample.fastq
	
	#Output directory
	output=$directory/ribose-seq/results/$index/$sample/Alignment

	#Create directory
    	mkdir -p $output

	#Output files
	finalBAM=$output/$sample.bam
	statistics=$output/$sample-Statistics.txt
	
#############################################################################################################################
	#Trim FASTQ files based on quality and Illumina adapter content
	#java -jar $path/trimmomatic-0.36.jar SE -phred33 $fastq $output/QCtrimmed.fastq \
	#ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:10 SLIDINGWINDOW:5:15 MINLEN:$MIN

	#Trim UMI from 5' ends of reads (add UMI into read name for de-duplication step)
	#umitools trim --end 5 $output/QCtrimmed.fastq $UMI | gzip -c > $output/UMItrimmed.fastq
	
	#Trim UMI from 5' ends of reads (append UMI to read name for subsequent de-duplication step)
	umi_tools extract -I $output/QCtrimmed.fastq -p $UMI -L $output/log.file -S $output/UMItrimmed.fastq
	
	#Reverse complement reads (rNMP=reverse complement of 5' base)
	cat $output/UMItrimmed.fastq | seqtk seq -r - > $output/reverseComplement.fastq
		
	#Align reads to reference using Bowtie2 and output statistics
	#bowtie -m 1 $index $output/reverseComplement.fastq -S $output/temp.sam 2> $statistics
	#bowtie2 -x $index -U $output/reverseComplement.fastq -S $output/temp.sam 2> $statistics
	
	#Directly convert SAM file to sorted BAM file and create index for BAM file
	#samtools view -bS $output/temp.sam | samtools sort - -o $output/temp.bam && samtools index $output/temp.bam
	
	#Save only mapped reads to another sorted BAM file and create index for BAM file
	#samtools view -bF4 $output/temp.bam | samtools sort - -o $output/mapped.bam; samtools index $output/mapped.bam

	#Save only unmapped reads to another sorted BAM file and create index for BAM file
	#samtools view -bf4 $output/temp.bam | samtools sort - -o $output/unmapped.bam; samtools index $output/unmapped.bam
	
	#De-duplicate reads based on UMI and read position and create index for BAM file
	#umi_tools dedup -I $output/mapped.bam -S $finalBAM -L dedup.log && samtools index $finalBAM
	
	#De-duplicate reads based on UMI and create index for BAM file
	#umitools rmdup $output/mapped.bam $finalBAM > $output/temp.bed.gz && samtools index $finalBAM

	#Remove temporary files
	#rm -f $output/reverseComplement.fastq $output/temp.bam $output/temp.bam.bai temp.sam $output/mapped.bam
		
	#Notify user that alignment step is complete for which samples
	#echo "Alignment of $sample to $index reference genome is complete"
done
