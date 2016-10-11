#! /usr/bin/env bash

#Author: Alli Gombolay
#This script calculates nucleotide frequencies in sacCer2 genome (chr I-XVI, chr M, and 2micron).
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program (4_Background-Frequencies.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-r] 'reference genome' [-d] 'Ribose-seq directory' [-h]
          -r Reference genome of interest (i.e., sacCer2)
          -d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-r], [-d], and [-h])
while getopts "r:d:h" opt;
do
    case $opt in
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

#INPUT
#Location of file containing sequences of reference genome
FASTA=$directory/ribose-seq/scripts/$reference.fa

#OUTPUT
#Location of output "ribose-seq" backgroundNucleotideFrequencies directory
output="$directory/ribose-seq/results/Background-Nucleotide-Frequencies"

#Create output directory if it does not already exist
if [[ ! -d $output ]];
then
  mkdir $output
fi

#CALCULATION OF BACKGROUND NUCLEOTIDE FREQUENCIES

#Calculate frequencies of nucleotides in reference genome
genome="$output/$reference.genome.nucleotide.frequencies.tab"
python2.7 4_Background-Frequencies.py $FASTA --region-size-minimum 1 --region-size-maximum 3 --verbose > $genome

#Calculate frequencies of nucleotides in chrM of sacCer2 genome
chrM="$output/$reference.chrM.nucleotide.frequencies.tab"
python2.7 4_Background-Frequencies.py $FASTA --region-size-minimum 1 --region-size-maximum 3 --only-chrom chrM --verbose > $chrM

#Calculate frequencies of nucleotides in 2micron of sacCer2 genome 
plasmid="$output/$reference.2micron.nucleotide.frequencies.tab"
python2.7 4_Background-Frequencies.py $FASTA --region-size-minimum 1 --region-size-maximum 3 --only-chrom 2micron --verbose > $plasmid
