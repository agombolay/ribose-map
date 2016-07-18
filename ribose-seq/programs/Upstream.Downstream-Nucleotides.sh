#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: July 11, 2016
#This program determines 0-based coordinates of +/100 downstream/upstream nucleotides from each rNMP

#COMMAND LINE OPTIONS

#Name of the program (Nucleotide-Coordinates.sh)
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

#Obtain 0-based coordinates of +/- 100 downstream/upstream nucleotides from each rNMP
for samples in ${bam[@]};
do

	#Extract sample names from filepaths
        filename=$(basename "${bam}")
        samples="${filename%.*}"

        #Extract input directories from filepaths
        inputDirectory=$(dirname "${bam}")

        #INPUT
	#Location of input BAM files
        input=$inputDirectory/$samples.bam

        #OUTPUT
        #Location of output "ribose-seq" alignment directory
        output=$directory/ribose-seq/results/$reference/$samples/alignment

	#Obtain positions of rNMPs (3’ end of each mapped read)
	bedtools genomecov -3 -bg -ibam $input > $samples.rNMPs.bothStrands.bed

	#Remove column containing coverage values
	awk '!($4="")' $samples.rNMPs.bothStrands.bed > temporary.bed

	#Change file back to original name
	mv temporary.bed $samples.rNMPs.bothStrands.bed

	#Make columns of BED file tab-delimited
	sed 's/ \+/\t/g' $samples.rNMPs.bothStrands.bed

	#Obtain coordinates of sacCer2 sequences that are +/- 100 bp downstream/upstream of each rNMP position:
	bedtools flank -i $samples.rNMPs.bothStrands.bed -g $directory/ribose-seq/reference/$reference.bed -b 100 > $samples.flanking.intervals.bed

	#Obtain sequences of sacCer2 coordinates from above that are +/- 100 bp downstream/upstream of each rNMP position:
	bedtools getfasta -fi $directory/ribose-seq/reference/$refernce.fa -bed $samples.flanking.intervals.bed -fo $samples.flanking.intervals.fasta

done
