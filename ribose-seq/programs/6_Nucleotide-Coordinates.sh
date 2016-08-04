#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: July 11, 2016
#This program determines 0-based coordinates of +/100 downstream/upstream nucleotides from each rNMP

#COMMAND LINE OPTIONS

#Name of the program (Nucleotide-Coordinates.sh)
program=$0

#Usage statement of the program
function usage () {
	echo "Usage: Nucleotide-Coordinates.sh [-i] 'BAM' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Filepaths of input BAM files ('/path/to/FS1.bam' etc.)
	-r One reference genome of interest ('sacCer2,' 'eColi,' 'hg38,' etc.)
	-d Directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
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
        output=$directory/ribose-seq/results/$reference/$samples/Nucleotide-Frequencies/Nucleotides

	#Create directory for output if it does not already exist
	if [[ ! -d $output ]];
	then
    		mkdir -p $output
	fi	

	#Location of output files
	positionsBoth=$output/$samples.rNMPs.bothStrands.bed

	flankingUpstreamIntervals=$output/$samples.flanking.upstream.intervals.bed
	flankingDownstreamIntervals=$output/$samples.flanking.downstream.intervals.bed

	flankingUpstreamSequences=$output/$samples.flanking.upstream.sequences.tab
        flankingDownstreamSequences=$output/$samples.flanking.downstream.sequences.tab

	temporary=$output/temporary.bed
	temporary2=$output/temporary2.bed

	#Obtain positions of rNMPs (3’ end of each mapped read)
	bedtools genomecov -3 -bg -ibam $input > $positionsBoth

	#Remove column containing coverage values
	awk '!($4="")' $positionsBoth > $temporary

	#Change file back to its original name
	mv $temporary $positionsBoth

	#Make columns of BED file tab-delimited
	sed 's/ \+/\t/g' $positionsBoth > $temporary2

	#Change file back to its original name
        mv $temporary2 $positionsBoth

	#Obtain coordinates of sacCer2 sequences that are 100 bp upstream of each rNMP position:
	bedtools flank -i $positionsBoth -g $directory/ribose-seq/reference/$reference.bed -l 100 -r 0 > $flankingUpstreamIntervals

	#Obtain coordinates of sacCer2 sequences that are 100 bp downstream of each rNMP position:
        bedtools flank -i $positionsBoth -g $directory/ribose-seq/reference/$reference.bed -l 0 -r 100 > $flankingDownstreamIntervals

	#Obtain sequences of sacCer2 coordinates from above that are 100 bp upstream of each rNMP position:
	bedtools getfasta -fi $directory/ribose-seq/reference/$reference.fa -bed $flankingUpstreamIntervals -tab -fo $flankingUpstreamSequences

	#Obtain sequences of sacCer2 coordinates from above that are 100 bp downstream of each rNMP position:
        bedtools getfasta -fi $directory/ribose-seq/reference/$reference.fa -bed $flankingDownstreamIntervals -tab -fo $flankingDownstreamSequences

done
