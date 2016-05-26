#! /usr/bin/env bash

output=$directory/ribose-seq/results/peaksSummary

if [[ ! -d $output ]]; then
    mkdir $output
fi

strands=("pos" "neg")

for strand in ${strands[@]}; do
    for align_mode in ${ALIGN_MODES[@]}; do

        samples=""
        names=""
        for sample in ${SAMPLES[@]}; do
            peakdir=$directory/ribose-seq/results/$samples/peaks
            exp_name="$sample.align.$align_mode.strand.$strand"
            peakbed="$peakdir/$exp_name""_peaks.narrowPeak"
            samples="$samples $peakbed"
            names="$names $sample"
        done

        result="$output/peaks_summary_table.align.$align_mode.strand.$strand.tab"

        bedtools multiinter -i $samples -names $names \
            -g $CHROM_SIZES \
            > $output
    done
done
