#! /usr/bin/env bash

#Author: Alli Gombolay
#Date: July 18, 2016
#This program outputs upstream and downstream flanking sequences into tabular format for processing

#COMMAND LINE OPTIONS

#Name of the program (Columns.sh)
program=$0

#Usage statement of the program
function usage () {
	echo "Usage: $program [-i] '/path/to/file1.tab etc.' [-l] 'location' [-r] 'reference' [-d] 'Ribose-seq directory' [-h]
	-i Filepaths of input tab files
	-l Upstream or Downstream location
	-r Reference genome of interest (i.e., sacCer2)
	-d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-i], [-l], [-r], [-d], and [-h])
while getopts "i:s:l:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) tab=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
        l ) location=$OPTARG ;;
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

for samples in ${tab[@]};
do

	#Extract sample names from filepaths
        filename=$(basename "${tab}")
        samples="${filename%.flanking.*}"

        #Extract input directories from filepaths
        inputDirectory=$(dirname "${tab}")

        #INPUT
	#Location of input tab files
	input=$inputDirectory/$samples.flanking.$location.sequences.tab

	#OUTPUT
	#Location of output "ribose-seq" Columns directory
	output=$directory/ribose-seq/results/$reference/$samples/nucleotideFrequencies/Nucleotides/Columns/$subset/$location

	if [[ ! -d $output ]]; then
    		mkdir -p $output
	fi

	#Location of output files
	sequences=$output/$samples.trimmed.$location.sequences.$subset.raw.txt
	columns=$output/$samples.trimmed.$location.sequences.$subset.columns.txt

	if [ $subset == "sacCer2" ];
	then
		cat $input > $sacCer2
	elif [ $subset == "chrM" ]
		grep 'chrM' $input > $chrM
	elif [ $subset == "nuclear" ]
		grep -v 'chrM' $input > $nuclear
	fi

	#Print sequences to new file
	awk -v "OFS=\t" '{print $2}' $genome > $sequences

	#Insert tabs between each nucleotide
	cat $sequences1 | sed 's/.../& /2g;s/./& /g' > $columns

		for i in {1..100};
		do
			awk -v field=$i '{ print $field }' $columns > $output/column.$i.$location.txt
		done
done

