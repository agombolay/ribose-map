#! /usr/bin/env bash

#Author: Alli Gombolay
#Date: July 18, 2016
#This program outputs upstream and downstream flanking sequences into tabular format for processing

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 7_Columns.sh [-i] 'TAB' [-s] 'Subset' [-l] 'Location' [-r] 'Reference' [-d] 'Directory' [-h]
	-i '/path/to/Nucleotide-Frequencies/Nucleotides/FS1.flanking.upstream.sequences.tab'
	-s Subset of reference genome of interest (sacCer2, hg38, eColi, nuclear, chrM, etc.)
	-l Location (Upstream or Downstream) specified in name of input TAB file shown above
	-r Name of reference genome folder in which to store output files ('sacCer2', etc.)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
}

#Use getopts function to create the command-line options ([-i], [-l], [-r], [-d], and [-h])
while getopts "i:s:r:d:h" opt;
do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        #i ) tab=($OPTARG) ;;
        i ) sample=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
        #l ) location=$OPTARG ;;
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

locations="upstream downstream"
tab="$sample.flanking.$location.sequences.tab"

for file in ${tab[@]};
do

	for location in ${locations[@]};
	do
		#Extract sample names from filepaths
        	#filename=$(basename "${tab}")
        	#samples="${filename%.flanking.*}"

        	#Extract input directories from filepaths
        	#inputDirectory=$(dirname "${tab}")

        	#INPUT
		#Location of input tab files
		#input=$inputDirectory/$samples.flanking.$location.sequences.tab

		file=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Nucleotides/$tab
	
		#OUTPUT
		#Location of output "ribose-seq" Columns directory
		output=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Nucleotides/Columns/$subset/$location

		if [[ ! -d $output ]]; then
    			mkdir -p $output
		fi

		if [[ ! -d $output/sequences ]]; then
                	mkdir -p $output/sequences
        	fi

		#Location of output files
		selection=$output/sequences/$sample.trimmed.$location.sequences.$subset.txt
		sequences=$output/sequences/$sample.trimmed.$location.sequences.$subset.raw.txt
		columns=$output/sequences/$sample.trimmed.$location.sequences.$subset.columns.txt

		if [ $subset == "sacCer2" ];
		then
			cat $input > $selection
		elif [ $subset == "chrM" ];
		then
			grep 'chrM' $input > $selection
		elif [ $subset == "nuclear" ];
		then
			grep -v 'chrM' $input > $selection
		fi

		#Print sequences to new file
		awk -v "OFS=\t" '{print $2}' $selection > $sequences

		#Insert tabs between each nucleotide
		cat $sequences | sed 's/.../& /2g;s/./& /g' > $columns

			for i in {1..100};
			do
				awk -v field=$i '{ print $field }' $columns > $output/$sample.column.$i.$location.$subset.txt
			done
	done
done
