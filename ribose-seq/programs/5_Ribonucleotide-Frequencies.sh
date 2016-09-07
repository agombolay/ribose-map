#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates the ribonucleotide frequencies located at 3' position of input BED file

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 5_Ribonucleotide-Frequencies.sh [-i] 'Sample' [-s] 'Subset' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of reference genome (sacCer2, hg38, eColi, nuclear, chrM, etc.)
	-r Name of reference genome folder in which to store output files (sacCer2, etc.)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
}

#Use getopts function to create the command-line options ([-i], [-s], [-r], [-d], and [-h])
while getopts "i:s:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) bed=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
	#r ) reference=$OPTARG ;;
	#d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ];
then
        exit
fi

#CALCULATE 3' NUCLEOTIDE FREQUENCIES

#Print only ribonucleotides (3' end of read (end for + strand and start for - strand))
cat $bed | while read -r line;
do

	#if [[ $subset == "sacCer2" && $line == *"+"* ]];
	#then
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,length,1)}' > temporary1
	#	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
	#	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt				
	
	#elif [[ $subset == "sacCer2" && $line == *"-"* ]];
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,0,1);}' > temporary1
        #	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
        #	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt
	
	#elif [[ $subset == "nuclear" && $line != *"chrM"* && $line == *"+"* ]];
	#then
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,length,1)}' > temporary1
	#	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
	#	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt
	#fi

	#elif [[ $subset == "nuclear" && $line != *"chrM"* && $line == *"-"* ]];
	#then
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,0,1);}' > temporary1
        #	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
        #	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt
	#fi

	#elif [[ $subset == "chrM" && $line == *"chrM"* && $line == *"+"* ]];
	#then
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,length,1)}' > temporary1
	#	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
	#	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt
	#fi

	#elif [[ $subset == "chrM" && $line == *"chrM"* && $line == *"-"* ]];
	#then
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,0,1);}' > temporary1
        #	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
        #	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt
	#fi

	#elif [[ $subset == "2micron" && $line == *"2micron"* && $line == *"+"* ]];
	#then
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,length,1)}' > temporary1
	#	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
	#	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt
	#fi

	#elif [[ $subset == "2micron" && $line == *"2micron"* && $line == *"-"* ]];
	#then
	#	awk -v "OFS=\t" '{print $5}' FS15.coordinates.bed | awk '{print substr($0,0,1);}' > temporary1
        #	awk -v "OFS=\t" '{print $6}' FS15.coordinates.bed > temporary2
        #	paste temporary1 temporary2 > $directory/$sample.List.$subset.txt
	#fi

#file="alli.txt"

#for file in $file;
#do
#		A=$(grep -o 'A' $file | wc -l)
#		C=$(grep -o 'C' $file | wc -l)
#		G=$(grep -o 'G' $file | wc -l)
#		T=$(grep -o 'T' $file | wc -l)

#		total=$(($A+$C+$G+$T))
	
#		A_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $A/$total`")
#		C_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $C/$total`")
#		G_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $G/$total`")
#		T_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $T/$total`")
		
#		A_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $A_frequency/$A_background_frequency`")
#        	C_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $C_frequency/$C_background_frequency`")
#        	G_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $G_frequency/$G_background_frequency`")
#        	T_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $T_frequency/$T_background_frequency`")

done
