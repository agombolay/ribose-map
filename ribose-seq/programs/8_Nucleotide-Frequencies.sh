#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates the frequencies of +/- 100 downstream/upstream nucleotides from ribonucleotides 

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 8_Nucleotide-Frequencies.sh [-i] 'Sample' [-s] 'Subset' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of reference genome (sacCer2, hg38, eColi, nuclear, chrM, etc.)
	-r Name of reference genome folder in which to store output files (sacCer2, etc.)
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

#Remove old .txt files
rm $output2/*.txt

locations="upstream downstream"

for location in ${locations[@]};
do

	A_normalized_frequencies=$output2/A_normalized_frequencies.$subset.$location.txt
	C_normalized_frequencies=$output2/C_normalized_frequencies.$subset.$location.txt
	G_normalized_frequencies=$output2/G_normalized_frequencies.$subset.$location.txt
	T_normalized_frequencies=$output2/T_normalized_frequencies.$subset.$location.txt
	Normalized_Frequencies=$output2/$sample.Normalized_Frequencies.$subset.$location.txt
		
	input=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Nucleotides/Columns/$subset/$location/$sample*.txt
	
	for file in ${input[@]};
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

		echo $A_normalized_frequency >> $A_normalized_frequencies
		echo $C_normalized_frequency >> $C_normalized_frequencies
		echo $G_normalized_frequency >> $G_normalized_frequencies
		echo $T_normalized_frequency >> $T_normalized_frequencies

		if [ -e "$Normalized_Frequencies" ]; then
    			rm $Normalized_Frequencies
		fi

		paste $A_normalized_frequencies $C_normalized_frequencies $G_normalized_frequencies $T_normalized_frequencies >> $Normalized_Frequencies
	done
done

#Create new folder, "data," to store output nucleotide frequency data files
data=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/$sample-Data

#Create folder only if it does not already exist
if [[ ! -d $data ]];
then
	mkdir $data    
fi

#Copy data files to new folder
cp $Normalized_Frequencies $data
cp $directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Ribonucleotides/$sample.Ribonucleotide.Frequencies.*.txt $data
cp $directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Ribonucleotides/$sample.Ribonucleotides-List.*.txt $data
