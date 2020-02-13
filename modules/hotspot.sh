#!/usr/bin/env bash

#Load config file
. "$1"

#Create output directory and remove any old files
output=$repository/results/$sample/hotspot$quality
rm -r $output; mkdir -p $output

######################################################################################################################################################################################################

for region in $other "chromosomes"; do

	#Calculate index
	index=$(echo "$(wc -l < $repository/results/$sample/coordinate$quality/${sample}-$region.counts.tab)*0.99" | bc)

	#Test if index is integer or floating
	integer=$(echo $index | grep -E '[0-9]+\.[0]{2}' -)

######################################################################################################################################################################################################

	if [[ $integer ]]; then
		
		#Calculate line numbers
		line1=$(echo $index | awk '{print int($1)}')
		line2=$(echo "$line1+1" | bc)

		#Get counts of rNMPs for lines one and two
		value1=$(head -${line1} $output/${sample}-$region.counts.tab | tail -1 | awk '{ print $5 }')
		value2=$(head -${line2} $output/${sample}-$region.counts.tab | tail -1 | awk '{ print $5 }')

		#Calculate 99th percentile (average of counts)
		percentile=$(echo "scale=2; ($value1 + $value2)/2" | bc -l)

	else
		
		#Round up value of index to get line number
		line=$(echo $index | awk '{print int($1 + 1)}')

		#Calculate 99th percentile (count at the line)
		percentile=$(head -${line} $output/${sample}-$region.counts.tab | tail -1 | awk '{ print $5 }')
	
	fi

######################################################################################################################################################################################################

	#Save top 1% of rNMP coordinates
	awk -v "OFS=\t" -v "x=$percentile" '{if ($5 >= x) print $0}' $output/${sample}-$region.counts.tab | awk -v "OFS=\t" '{ print $1, $2, $3, ".", ".", $4 }' - > $output/${sample}-$region.top.bed

######################################################################################################################################################################################################

	for file in $(ls $output/${sample}-$region.top.bed); do

		#Get sequence of rNMPs and the 3 nucleotides upstream and downstream from them
		bedtools slop -s -i $file -g $output/$(basename $fasta .fa).bed -b 3 > $output/$(basename $file .bed).slop.bed
	
		bedtools getfasta -s -fi $fasta -bed $output/$(basename $file .bed).slop.bed | awk '/>/{$0 = ">" ++i substr($0, 2)} 1' - > $output/$(basename $file .bed).flank.txt

		sites=$(grep -c "^>" $output/$(basename $file .bed).flank.txt)

		meme $output/$(basename $file .bed).flank.txt -o $output/meme-$(basename $file .bed) -dna -minw 7 -nsites $sites -brief 1000000

	done
done

rm $output/*.{tab,txt}
