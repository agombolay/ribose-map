#! /usr/bin/env bash

#Author: Alli Gombolay
#This program calculates the frequencies of the nucleotides located in the assembled samples.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

#COMMAND LINE OPTIONS

#Name of the program (6_Experimental-Frequencies.sh)
program=$0

#Usage statement of the program
function usage () {
        echo "Usage: $program [-i] 'sample1 etc.' [-r] 'reference genome' [-d] 'Ribose-seq directory' [-h]
        -i Sample names of input BAM files (i.e, sample1 for sample1.bam)
        -r File containing sizes in base pairs of chromosomes (i.e, sacCer2)
        -d Location of user's local Ribose-seq directory"
}

#Use getopts function to create the command-line options ([-d] and [-h])
while getopts "i:r:d:h" opt;
do
    case $opt in
	 #Specify input as arrays to allow multiple input arguments
        i ) samples=($OPTARG) ;;
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

# Mononucleotides, dinucleotides, and trinucleotides
sizes="1 2 3"

modes=("all" "only-mitochondria" "no-mitochondria" "only-2micron")

arguments=("" "--only-chrom chrM" "--ignore-chrom chrM" "--only-chrom 2micron")

#""= Entire genome ("all")
#"--only-chrom chrM"= Only chrM ("only-mito")
#"--only-chrom 2micron"=Only 2micron plasmid ("only-2micron")
#"--ignore-chrom chrM"= Nuclear chromosomes and 2micron plasmid ("no-mito")

input=$directory/ribose-seq/results/$samples/alignment

output=$directory/ribose-seq/results/$samples/nucleotideFrequencies

if [[ ! -d $output ]]; then
    mkdir -p $output
fi

offset_minimum=-100
offset_maximum=100

BAM=$input/$samples.bam

#Location of file containing sequences of reference genome
FASTA=$directory/ribose-seq/reference/$reference.fa

for index in ${!modes[@]};
do

        mode=${modes[$index]}
        
        argument=${arguments[$index]}

        tables="$output/$samples.$mode.nucleotideFrequencies.tab"
        
        if [[ -f $tables ]];
        then
            rm -f $tables
        fi

        if [[ $mode == "only-mitochondria" ]];
        then
            BackgroundFrequencies="$directory/ribose-seq/results/backgroundNucleotideFrequencies/$reference.chrM.nucleotide.frequencies.tab"
        
        elif [[ $mode == "only-2micron" ]];
        then
            BackgroundFrequencies="$directory/ribose-seq/results/backgroundNucleotideFrequencies/$reference.2micron.nucleotide.frequencies.tab"
        
        else
            BackgroundFrequencies="$directory/ribose-seq/results/backgroundNucleotideFrequencies/$reference.genome.nucleotide.frequencies.tab"
        fi

        #Signals need to be reverse complemented since sequence is reverse complement of the captured strand
        for size in $sizes;
        do
<<<<<<< HEAD
        	/projects/home/agombolay3/data/bin/Python-2.7.11/python 5_Experimental-Frequencies.py $BAM $FASTA \
		--region-size $size $arguments --revcomp-strand --background-freq-table $BackgroundFrequencies \
		--offset-min $offset_minimum --offset-max $offset_maximum --verbose >> $tables
=======
            /projects/home/agombolay3/data/bin/Python-2.7.11/python 5_Experimental-Frequencies.py $BAM $FASTA
            --region-size $size $arguments --revcomp-strand --background-freq-table $BackgroundFrequencies
            --offset-min $offset_minimum --offset-max $offset_maximum --verbose >> $tables
>>>>>>> 00e51bf587afd72b8c69ce7f5fd7a4a8d652fad0
        done
    
done
