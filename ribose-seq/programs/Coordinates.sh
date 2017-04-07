#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program determines the coordinates of rNMPs (3' position of aligned reads)

#Usage statement
function usage () {
	echo "Usage: Coordinates.sh [-i] 'Sample(s)' [-s] 'Subset' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Sample name(s) (FS1, FS2, FS3 etc.)
	-s Subset of genome (genome, nuclear, mito)
	-r Reference genome (sacCer2, pombe, ecoli, mm9, hg38, etc.)
	-d Directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Command-line options
while getopts "i:s:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Allow only one input argument
	s ) subset=$OPTARG ;;
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Determine coordinates
for sample in ${sample[@]}; do

#############################################################################################################################
	#Input file
	bam=$directory/ribose-seq/results/$reference/$sample/Alignment/$sample.bam
	
	#Output directory
	output=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset

	#Create directory
	mkdir -p $output
	
	#Remove older versions of files
	rm -f $output/{*.txt,*.bed,*.fa,*.fq}
	
	#Output files
	reads=$output/$sample.read-information.$subset.txt
	coordinates=$output/$sample.rNMP-coordinates.$subset.bed

#############################################################################################################################
	#STEP 1: Extract sequences from BAM alignment file

	#Convert BAM to FASTA file then extract sequences from FASTA
	samtools bam2fq $bam | seqtk seq -A - | grep -v '>' - > temporary1

#############################################################################################################################
	#STEP 2: Obtain rNMP coordinates from aligned reads

	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > temporary2
	
	#Extract read coordinates, sequences, and strand information
	paste temporary2 temporary1 | awk -v "OFS=\t" '{print $1, $2, $3, $4, $6, $7}' > $reads
	
	#Obtain coordinates of rNMPs located on positive strand of DNA
	positiveReads=$(awk -v "OFS=\t" '$5 == "+" {print $1, ($3 - 1), $3, " ", " ", $5}' $reads)
	
	#Obtain coordinates of rNMPs located on negative strand of DNA
	negativeReads=$(awk -v "OFS=\t" '$5 == "-" {print $1, $2, ($2 + 1), " ", " ", $5}' $reads)
		
	#Filter coordinates by region
	if [ $subset == "nuclear" ]; then
		#Combine +/- coordinates and select only nuclear DNA
		cat <(echo "$positiveReads") <(echo "$negativeReads") \
		| grep -v -E '(chrM|MT|AB325691|chrEBV|chrUN*|*random)' - > temporary3
	
	elif [ $subset == "mito" ]; then
		#Combine +/- coordinates and select only mitochondrial DNA
		cat <(echo "$positiveReads") <(echo "$negativeReads") | grep -E '(chrM|MT)' - > temporary3
	
	else
		#Combine all +/- rNMP coordinates
		cat <(echo "$positiveReads") <(echo "$negativeReads") > temporary3
	fi
	
	#Sort ribonucleotide coordinates
	sort -k1,1 -k2,2n temporary3 > $coordinates
done

#Remove temporary files
rm temporary1 temporary2 temporary3
