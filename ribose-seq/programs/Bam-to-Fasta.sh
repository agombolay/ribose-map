#!/usr/bin/env bash

#Author: Alli Gombolay
#This program removes UMI's from reads, aligns reads to reference genome, and de-duplicates reads
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program (BAM-to-FASTA.sh)
program=$0

#Usage statement of the program
function usage () {
	echo "Usage: $program [-i] '/path/to/file1.bam etc.' [-r] 'reference' [-d] 'Ribose-seq directory' [-h]
	-i Filepaths of input BAM files
	-r Reference genome of interest (i.e., sacCer2)
	-d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-i], [-r], [-d], and [-h])
while getopts "i:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) bam=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

for samples in ${bam[@]};
do

	#Extract sample names from filepaths
	filename=$(basename "${samples}")
	samples="${filename%.*}"

	#Extract input directories from filepaths
	inputDirectory=$(dirname "${bam}")

	#Location of input BAM files
	input=$inputDirectory/$samples.bam

	#OUTPUT
	#Location of output "ribose-seq" nucleotideFrequencies directory
	output=$directory/ribose-seq/results/$reference/$samples/Nucleotide-Frequencies/Ribonucleotides

	fastq=$output/$samples.fastq
	fasta=$output/$samples.fasta

	samtools bam2fq $input > $fastq
	seqtk seq -A $fastq > $fasta

done
