#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates the frequencies of +/- 100 downstream/upstream nucleotides from ribonucleotides 

#COMMAND LINE OPTIONS

#Name of the program (Nucleotide-Frequencies.sh)
program=$0

#Usage statement of the program
function usage () {
	echo "Usage: $program [-i] 'TXT' [-s] 'Subset' [-l] 'Location' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Filepaths of text files of each column of data ('/path/to/*.txt' etc.)
	-s Genetic subset of genome of interest ('sacCer2,' 'nuclear,' 'chrM,' etc.)
	-l Location (Upstream or Downstream) specified in name of input TXT files above
	-r Name of reference genome folder in which to store output files (i.e., sacCer2)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
}

#Use getopts function to create the command-line options ([-i], [-s], [-l], [-n], [-r], [-d], and [-h])
while getopts "i:s:l:n:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) txt=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
	l ) location=$OPTARG ;;
	n ) sample=$OPTARG ;;
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

#Extract sample name from filepath
filename=$(basename "${txt}")
sample="${filename%.column.*}"

#INPUT
#Location of FASTA file
fasta=$directory/ribose-seq/reference/$subset.fa

#OUTPUT
#Location of output directory
output1=$directory/ribose-seq/results/Background-Nucleotide-Frequencies
output2=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Nucleotides

for file in $fasta;
do
  	A=$(grep -v '>' $file | grep -o 'A' - | wc -l)
        C=$(grep -v '>' $file | grep -o 'C' - | wc -l)
        G=$(grep -v '>' $file | grep -o 'G' - | wc -l)
        T=$(grep -v '>' $file | grep -o 'T' - | wc -l)

        total=$(($A+$C+$G+$T))

        A_background_frequency=$(bc <<< "scale = 4; `expr $A/$total`")
        C_background_frequency=$(bc <<< "scale = 4; `expr $C/$total`")
        G_background_frequency=$(bc <<< "scale = 4; `expr $G/$total`")
        T_background_frequency=$(bc <<< "scale = 4; `expr $T/$total`")
	
	echo "A Background Frequency: $A_background_frequency" >> $output1/Background-Frequencies.txt
	echo "C Background Frequency: $C_background_frequency" >> $output1/Background-Frequencies.txt
	echo "G Background Frequency: $G_background_frequency" >> $output1/Background-Frequencies.txt
	echo "T Background Frequency: $T_background_frequency" >> $output1/Background-Frequencies.txt

done

for file in ${txt[@]};
do
	A=$(grep -v '>' $file | grep -o 'A' - | wc -l)
	C=$(grep -v '>' $file | grep -o 'C' - | wc -l)
	G=$(grep -v '>' $file | grep -o 'G' - | wc -l)
	T=$(grep -v '>' $file | grep -o 'T' - | wc -l)

	total=$(($A+$C+$G+$T))
	
	A_frequency=$(bc <<< "scale = 4; `expr $A/$total`")
	C_frequency=$(bc <<< "scale = 4; `expr $C/$total`")
	G_frequency=$(bc <<< "scale = 4; `expr $G/$total`")
	T_frequency=$(bc <<< "scale = 4; `expr $T/$total`")

	A_normalized_frequency=$(bc <<< "scale = 4; `expr $A_frequency/$A_background_frequency`")
        C_normalized_frequency=$(bc <<< "scale = 4; `expr $C_frequency/$C_background_frequency`")
        G_normalized_frequency=$(bc <<< "scale = 4; `expr $G_frequency/$G_background_frequency`")
        T_normalized_frequency=$(bc <<< "scale = 4; `expr $T_frequency/$T_background_frequency`")

	echo $A_normalized_frequency >> $output2/A_normalized_frequencies.$subset.$location.txt
	echo $C_normalized_frequency >> $output2/C_normalized_frequencies.$subset.$location.txt
	echo $G_normalized_frequency >> $output2/G_normalized_frequencies.$subset.$location.txt
	echo $T_normalized_frequency >> $output2/T_normalized_frequencies.$subset.$location.txt

	paste $output2/A_normalized_frequencies.$subset.$location.txt $output2/C_normalized_frequencies.$subset.$location.txt $output2/G_normalized_frequencies.$subset.$location.txt $output2/T_normalized_frequencies.$subset.$location.txt > $output2/Normalized-Frequencies.$subset.$location.txt

done
