#!/usr/bin/env bash
#Author: Alli Gombolay
#This program creates a .txt file containing the output nucleotide frequencies data values for plotting

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

#Print values -100 to 100
seq -100 1 100 > temporary1.txt

#Combine upstream and downstream normalized nucleotide frequency files and ribonucleotide frequency files
cat $sample.Normalized_Frequencies.$subset.upstream.txt $sample.Ribonucleotide.Frequencies.$subset.txt $sample.Normalized_Frequencies.$subset.downstream.txt >> temporary2.txt

output=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/$sample-Data/$sample.Nucleotide_Frequencies.$subset.txt

#Merge two files into final .txt file
paste temporary1.txt temporary2.txt > temporary3.txt

#Add Header to beginning of .txt file 
echo "Position A C G U/T" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5}' | cat - temporary3.txt > temp && mv temp temporary3.txt

#Make sure data values are arranged in columns
column -t temporary3.txt > $output

rm temporary1.txt temporary2.txt temporary3.txt
