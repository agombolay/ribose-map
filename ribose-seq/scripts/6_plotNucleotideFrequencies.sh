ignore_modes=("all" "only-mito" "no-mito" "only-2micron")

offset_values=(100 50 15)

output=$directory/ribose-seq/results/$sample/plots
subPlotDirectory="$output/nucleotideFrequencies"

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
    
    sampleID="$sample.subset-$ignore_mode"
    
    Rscript $RSCRIPTS/nuc.freqs.R -n "$sampleID" -d $subPlotDirectory --offsetmax $offset_values $tables
  done

done
