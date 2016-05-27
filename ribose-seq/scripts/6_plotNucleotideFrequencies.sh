# for nuc freq plots
offset_maxs=(100 50 15)

subplotdir="$plotdir/nuc_freqs"

if [[ ! -d $subplotdir ]];
then
  mkdir -p $subplotdir
fi

for ig_idx in ${!ignore_modes[@]};
do
  ignore_mode=${ignore_modes[$ig_idx]}

    for offset_max in ${offset_maxs[@]};
    do
      counts="$results/nuc_freqs/$sample.align.$align_mode.ignore.$ignore_mode.nuc_freqs.tab"
      sampleid="$sample.align-$align_mode.subset-$ignore_mode"

      Rscript $RSCRIPTS/nuc.freqs.R -n "$sampleid" -d $subplotdir --offsetmax $offset_max $counts
    done
