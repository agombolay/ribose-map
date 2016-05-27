output=$directory/ribose-seq/data/$sample/plots

# for nuc freq plots
offset_maxs=(100 50 15)

subplotdir="$output/nucleotideFrequencies"

if [[ ! -d $subplotdir ]];
then
  mkdir -p $subplotdir
fi

for index in ${!ignore_modes[@]};
do

  ignore_mode=${ignore_modes[$index]}

  for offset_max in ${offset_maxs[@]};
  do
    tables="$directory/ribose-seq/results/nucleotideFrequencies/$sample.ignore.$ignore_mode.nucleotideFrequencies.tab"
    sampleid="$sample.align-$align_mode.subset-$ignore_mode"
    Rscript $RSCRIPTS/nuc.freqs.R -n "$sampleid" -d $subplotdir --offsetmax $offset_max $counts
  done

done
