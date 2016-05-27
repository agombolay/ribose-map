#! /usr/bin/env bash

#BSUB -J nuc.freqs[1-13]
#BSUB -e nuc.freqs.%J.%I.err
#BSUB -o nuc.freqs.%J.%I.out
#BSUB -q normal
#BSUB -P storici

<<DOC
Calculate nucleotide frequencies
DOC

set -o nounset -o pipefail -o errexit -x

source $CONFIG
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# mono, di and trinucleotides
sizes="1 2 3"

if [[ $ASSEMBLY == "sacCer2" ]]; then
    ignore_modes=("all" "only-mito" "no-mito" "only-2micron")
    ignore_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM"
                 "--only-chrom 2micron")
else
    ignore_modes=("all" "only-mito" "no-mito")
    ignore_args=("" "--only-chrom chrM"
                 "--ignore-chrom chrM")
fi

aligndir=$RESULT/$sample/alignment
results=$RESULT/$sample/nuc_freqs
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

# need to be in BIN to run module
cd $BIN

offset_min=-100
offset_max=100

for aln_idx in ${!ALIGN_MODES[@]}; do
    align_mode=${ALIGN_MODES[$aln_idx]}
    BAM=$aligndir/$sample.align.$align_mode.bam

    for ig_idx in ${!ignore_modes[@]}; do

        ignore_mode=${ignore_modes[$ig_idx]}
        ignore_arg=${ignore_args[$ig_idx]}

        output="$results/$sample.align.$align_mode.ignore.$ignore_mode.nuc_freqs.tab"
        if [[ -f $output ]]; then
            rm -f $output
        fi

        if [[ $ignore_mode == "only-mito" ]]; then
            BKGD_FREQS="$RESULT/background_nuc_freqs/$ASSEMBLY.chrM.nuc.freqs.tab"
        elif [[ $ignore_mode == "only-2micron" ]]; then
            BKGD_FREQS="$RESULT/background_nuc_freqs/$ASSEMBLY.2micron.nuc.freqs.tab"
        else
            BKGD_FREQS="$RESULT/background_nuc_freqs/$ASSEMBLY.genome.nuc.freqs.tab"
        fi

        # signals need to be reverse complemented because the sequence is
        # the reverse complement of the captured strand
        for size in $sizes; do
            python -m modmap.nuc_frequencies \
                $BAM $FASTA \
                --region-size $size \
                $ignore_arg \
                --revcomp-strand \
                --background-freq-table $BKGD_FREQS \
                --offset-min $offset_min \
                --offset-max $offset_max \
                --verbose \
                >> $output
        done
    done
done
