ignore_modes=("all" "only-mito" "no-mito" "only-2micron")

offset_values=(100 50 15)

output=$directory/ribose-seq/results/$samples/plots
subPlotDirectory="$output/nucleotideFrequencies"

input=$directory/ribose-seq/results/$samples/nucleotideFrequencies/

if [[ ! -d $subplotdir ]];
then
  mkdir -p $subplotdir
fi

for index in ${!ignore_modes[@]};
do

  ignore_mode=${ignore_modes[$index]}

  for offset_max in ${offset_maxs[@]};
  do
    sampleID="$sample.subset-$ignore_mode"
    tables="input/$sample.ignore.$ignore_mode.nucleotideFrequencies.tab"
    Rscript nucleotideFrequencies.R -n "$sampleID" -d $subPlotDirectory --offsetmax $offset_values $tables
  done

done
