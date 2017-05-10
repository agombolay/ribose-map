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
	unmappedBAM=$output/$sample-unmapped.bam; unmappedFASTQ=$output/$sample-unmapped.fastq;
	statistics=$output/$sample-Statistics.txt; BED=$output/$sample.bed.gz; finalBAM=$output/$sample.bam
	tempSAM=$output/$sample-temp.sam; tempBAM=$output/$sample-temp.bam; mappedBAM=$output/$sample-mapped.bam;
	extracted=$output/$sample.UMI-extracted.fastq.gz; reverseComplement=$output/$sample-reverseComplement.fastq;
	
#############################################################################################################################
	#Trim FASTQ files based on quality and Illumina adapter content
	#java -jar $path/trimmomatic-0.36.jar SE -phred33 $fastq $output/$sample-trimmed.fastq \
	#ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:10 SLIDINGWINDOW:5:15 MINLEN:$MIN

	#Trim UMI from 5' ends of reads (add UMI into read name for further processing)
	#umitools trim --end 5 $output/$sample-trimmed.fastq $UMI | gzip -c > $umiTrimmed

	#umi_tools extract -I $output/$sample-trimmed.fastq -p $UMI -L log.file -S $umiTrimmed 
	
	#Reverse complement reads
	#zcat $umiTrimmed | seqtk seq -r - > $reverseComplement
	
	#Align reads to reference genome using Bowtie2
	#bowtie -m 1 $index $reverseComplement -S $tempSAM 2> $statistics
	bowtie2 -x $index -U $reverseComplement -S $tempSAM 2> $statistics
	
	#Directly convert SAM file to sorted BAM file and create index for BAM file
	samtools view -bS $tempSAM | samtools sort - -o $tempBAM && samtools index $tempBAM
	
	#Convert SAM file to BAM and sort temp BAM file
	#-S: Input=SAM; -h: header; -u: Output=uncompressed BAM
	#samtools view -Shu $tempSAM | samtools sort - -o $tempBAM
	
	samtools view -bF4 $tempBAM | samtools sort - -o $mappedBAM; samtools index $mappedBAM

	#Save mapped reads to BAM
	#"-F4": Output only mapped reads
	#samtools view -bF4 $tempBAM > $mappedBAM
		
	#Save unmapped reads to FASTQ
	#"-f4": Output only unmapped reads
	#samtools view -bf4 $tempBAM > $unmappedBAM
	#bamToFastq -i $unmappedBAM -fq $unmappedFASTQ
	
	#De-duplicate reads by saving one per UMI
	#umitools rmdup $mappedBAM $finalBAM > $BED
	
	#De-duplicate reads based on UMI
	#umi_tools dedup -I $mapped2BAM -S $finalBAM -L dedup.log
	
	#Create index file
	#samtools index $finalBAM
		
	#Remove temporary files
	#rm -f $mappedBAM $tempBAM $tempBAM.bai $tempSAM
		
	#Notify user that alignment step is complete for which samples
	#echo "Alignment of $sample to $index reference genome is complete"
done
