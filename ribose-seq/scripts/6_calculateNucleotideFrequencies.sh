#! /usr/bin/env bash

#Author: Alli Gombolay
#This program calculates the frequencies of the nucleotides located in the assembled samples.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

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

BAM=$input/$sample.bam

#Location of file containing sequences of reference genome
FASTA=$directory/ribose-seq/data/reference/$reference.fa

for index in ${!modes[@]};
do

        mode=${modes[$index]}
        
        argument=${arguments[$index]}

        tables="$output/$sample.$mode.nucleotideFrequencies.tab"
        
        if [[ -f $tables ]];
        then
            rm -f $tables
        fi

        if [[ $mode == "only-mitochondria" ]];
        then
            BackgroundFrequencies="$output/backgroundNucleotideFrequencies/chrM.nucleotide.frequencies.tab"
        
        elif [[ $mode == "only-2micron" ]];
        then
            BackgroundFrequencies="$output/backgroundNucleotideFrequencies/2micron.nucleotide.frequencies.tab"
        
        else
            BackgroundFrequencies="$output/backgroundNucleotideFrequencies/genome.nucleotide.frequencies.tab"
        fi

        #Signals need to be reverse complemented since sequence is reverse complement of the captured strand
        for size in $sizes;
        do
            python2.7 6_calculateNucleotideFrequencies.py $BAM $FASTA --region-size $size $arguments --revcomp-strand \
            --background-freq-table $BackgroundFrequencies --offset-min $offset_minimum --offset-max $offset_maximum \
            --verbose >> $output
        done
    
done
