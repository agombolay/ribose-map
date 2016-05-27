#! /usr/bin/env bash

# mono, di and trinucleotides
sizes="1 2 3"


ignore_modes=("all" "only-mito" "no-mito" "only-2micron")
ignore_args=("" "--only-chrom chrM" "--ignore-chrom chrM" "--only-chrom 2micron")

input=$directory/ribose-seq/data/$sample/alignment

output=$directory/ribose-seq/data/$sample/nucleotideFrequencies

if [[ ! -d $output ]]; then
    mkdir -p $output
fi

offset_min=-100
offset_max=100


BAM=$input/$sample.bam

for index in ${!ignore_modes[@]}; do

        ignore_mode=${ignore_modes[$index]}
        
        ignore_arg=${ignore_args[$index]}

        tables="$output/$sample.ignore.$ignore_mode.nucleotideFrequencies.tab"
        
        if [[ -f $tables ]];
        then
            rm -f $tables
        fi

        if [[ $ignore_mode == "only-mito" ]];
        then
            BKGD_FREQS="$RESULT/background_nuc_freqs/chrM.nuc.freqs.tab"
        
        elif [[ $ignore_mode == "only-2micron" ]];
        then
            BKGD_FREQS="$RESULT/background_nuc_freqs/2micron.nuc.freqs.tab"
        
        else
            BKGD_FREQS="$RESULT/background_nuc_freqs/genome.nuc.freqs.tab"
        fi

        #Signals need to be reverse complemented since the sequence is the reverse complement of the captured strand
        for size in $sizes;
        do
            python -m modmap.nuc_frequencies $BAM $FASTA --region-size $size $ignore_arg --revcomp-strand \
            --background-freq-table $BKGD_FREQS --offset-min $offset_min --offset-max $offset_max --verbose \
            >> $output
        done
    
done
