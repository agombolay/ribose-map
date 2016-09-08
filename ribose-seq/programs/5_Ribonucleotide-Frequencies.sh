#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates the ribonucleotide frequencies located at 3' position of input BED file

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 5_Ribonucleotide-Frequencies.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of genome (sacCer2, nuclear, chrM)
	-r Reference genome assembly (sacCer2, etc.)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
}

#Use getopts function to create the command-line options ([-i], [-s], [-r], [-d], and [-h])
while getopts "i:s:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
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

bed=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Ribonucleotides/$sample.coordinates.bed

#CALCULATE 3' NUCLEOTIDE FREQUENCIES
#Print only ribonucleotides of genome subset to output file

#Whole genome subset
if [[ $subset == "sacCer2" ]];
then
	cut -d'	' -f4,5 $bed > List.$subset.temp
#Nuclear subset
elif [[ $subset == "nuclear" ]];
then
	grep -v 'chrM' $bed | cut -d'	' -f4,5 - > List.$subset.temp
#Mitochondria subset
elif [[ $subset == "mitochondria" ]];
then
	grep 'chrM' $bed | cut -d'	' -f4,5 - > List.$subset.temp
fi

#Print only ribonucleotides (3' end of read (end for + strand and start for - strand)) to output file

#Print ribonucleotides for positive strands (located at end of sequence)
awk '$2 == "+" { print substr( $0, length($0) - 2, length($0) ) }' List.$subset.temp > List.$subset.positive.temp

#Print ribonucleotides for negative strands (located at start of sequence)
awk '$2 == "-" {print substr($0,0,1), $2;}' List.$subset.temp > List.$subset.negative.temp

#Combine output files generated from above into one final output file
paste List.$subset.positive.temp List.$subset.negative.temp > List.$subset.txt

A_ribonucleotide_count=$(grep -o 'A' List.$subset.txt | wc -l)
C_ribonucleotide_count=$(grep -o 'C' List.$subset.txt | wc -l)
G_ribonucleotide_count=$(grep -o 'G' List.$subset.txt | wc -l)
U_ribonucleotide_count=$(grep -o 'T' List.$subset.txt | wc -l)

total=$(($A_ribonucleotide_count+$C_ribonucleotide_count+$G_ribonucleotide_count+$U_ribonucleotide_count))
	
A_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $A_ribonucleotide_count/$total`")
C_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $C_ribonucleotide_count/$total`")
G_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $G_ribonucleotide_count/$total`")
U_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $U_ribonucleotide_count/$total`")
		
#A_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $A_ribonucleotide_frequency/$A_background_frequency`")
#C_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $C_ribonucleotide_frequency/$C_background_frequency`")
#G_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $G_ribonucleotide_frequency/$G_background_frequency`")
#U_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $U_ribonucleotide_frequency/$T_background_frequency`")

#echo $A_ribonucleotide_normalized_frequency
#echo $C_ribonucleotide_normalized_frequency
#echo $G_ribonucleotide_normalized_frequency
#echo $U_ribonucleotide_normalized_frequency

echo $A_ribonucleotide_frequency
