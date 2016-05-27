#! /usr/bin/env bash

#Author: Alli Gombolay
#This program calculates the frequencies of the nucleotides located in the assembled samples.
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

# mononucleotides, dinucleotides, and trinucleotides
sizes="1 2 3"

ignore_modes=("all" "only-mito" "no-mito" "only-2micron")

ignore_arguments=("" "--only-chrom chrM" "--ignore-chrom chrM" "--only-chrom 2micron")

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

for index in ${!ignore_modes[@]};
do

        ignore_mode=${ignore_modes[$index]}
        
        ignore_arguments=${ignore_args[$index]}

        tables="$output/$sample.ignore.$ignore_mode.nucleotideFrequencies.tab"
        
        if [[ -f $tables ]];
        then
            rm -f $tables
        fi

        if [[ $ignore_mode == "only-mito" ]];
        then
            BKGD_FREQS="$RESULT/backgroundNucleotideFrequencies/chrM.nucleotide.frequencies.tab"
        
        elif [[ $ignore_mode == "only-2micron" ]];
        then
            BackgroundFrequencies="$RESULT/backgroundNucleotideFrequencies/2micron.nucleotide.frequencies.tab"
        
        else
            BackgroundFrequencies="$RESULT/backgroundNucleotideFrequencies/genome.nucleotide.frequencies.tab"
        fi

        #Signals need to be reverse complemented since sequence is reverse complement of the captured strand
        for size in $sizes;
        do
            python -m modmap.nuc_frequencies $BAM $FASTA --region-size $size $ignore_arguments --revcomp-strand \
            --background-freq-table $BackgroundFrequencies --offset-min $offset_minimum --offset-max $offset_maximum
            --verbose >> $output
        done
    
done
