#! /usr/bin/env bash

#Author: Alli Gombolay
#Date: August 2, 2016
#This program allows the user to view align data from input BAM file using "samtools tview" program

#COMMAND LINE OPTIONS

#Name of the program (View-Alignment.sh)
program=$0

#Usage statement of the program
function usage () {
	echo "Usage: $program [-i] '/path/to/file1.bam etc.' [-c] 'chromosome' [-p] 'nucleotide position' [-r] '/path/to/reference.fa' [-h]
	-i Input (i.e., /projects/home/agombolay3/data/repository/Ribose-seq-Project/ribose-seq/results/hg38/FS56/Alignment/FS56.bam)
	-c Chromosome of interest (i.e., chr1)
	-p Nucleotide position of interest (i.e., 95)
	-r Reference genome sequence file in FASTA format"
}

#Use getopts function to create the command-line options ([-i], [-c], [-p], [-r], and [-h])
while getopts "i:c:p:r:h" opt;
do
    case $opt in
	#Specify input as variable to allow only one input argument
	i ) bam=$OPTARG ;;
	c ) chromosome=$OPTARG ;;
	p ) position=$OPTARG ;;
	r ) reference=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

samtools tview -p $chromosome:$position $bam --reference $reference
