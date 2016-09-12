#!/usr/bin/env bash

#Author: Alli Gombolay
#This program converts input BAM files to FASTA files for viewing

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 3_Bam-to-Fasta.sh [-i] 'BAM' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Filepaths of BAM files ('/path/to/FS1.bam' etc.)
	-r Name of reference genome folder in which to store output files ('sacCer2', etc.)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
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

	#Create directory for output if it does not already exist
	if [[ ! -d $output ]];
	then
    		mkdir -p $output
	fi
	
	fastq=$output/$samples.fastq
	fasta=$output/$samples.fasta

	samtools bam2fq $input > $fastq
	seqtk seq -A $fastq > $fasta

done
