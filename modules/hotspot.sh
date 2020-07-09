#!/usr/bin/env bash

#Load config file
. "$1"

#Create output directory and remove any old files
output=$repository/results/$sample/hotspot$quality
rm -r $output; mkdir -p $output

#############################################################################################################################################################################################################################################

for region in $other "chromosomes"; do
	
	for nuc in "A" "C" "G" "U" "Combined"; do
	
		#Calculate index
		index=$(echo "$(wc -l < $repository/results/$sample/coordinate$quality/${sample}-$region.$nuc.tab)*$percentile" | bc)

		#Test if index is integer or floating
		integer=$(echo $index | grep -E '[0-9]+\.[0]{2}' -)

#############################################################################################################################################################################################################################################

		if [[ $integer ]]; then
		
			#Calculate line numbers
			line1=$(echo $index | awk '{print int($1)}')
			line2=$(echo "$line1+1" | bc)

			#Get counts of rNMPs for lines one and two
			value1=$(head -${line1} $repository/results/$sample/coordinate$quality/${sample}-$region.$nuc.tab | tail -1 | awk '{ print $7 }')
			value2=$(head -${line2} $repository/results/$sample/coordinate$quality/${sample}-$region.$nuc.tab | tail -1 | awk '{ print $7 }')

			#Calculate Xth percentile (average of counts)
			threshold=$(echo "scale=2; ($value1 + $value2)/2" | bc -l)

		else
		
			#Round up value of index to get line number
			line=$(echo $index | awk '{print int($1 + 1)}')

			#Calculate Xth percentile (count at the line)
			threshold=$(head -${line} $repository/results/$sample/coordinate$quality/${sample}-$region.$nuc.tab | tail -1 | awk '{ print $7 }')
	
		fi
		
#############################################################################################################################################################################################################################################

		#Save top 1% of rNMP coordinates
		awk -v "OFS=\t" -v "x=$threshold" '{if ($7 >= x) print $0}' $repository/results/$sample/coordinate$quality/${sample}-$region.$nuc.tab > $output/${sample}-$region.$nuc.top.tab

#############################################################################################################################################################################################################################################

		#for file in $(ls $output/${sample}-$region.$nuc.top.tab); do

		#Get genomic coordinates of rNMPs and the 3 nucleotides up/downstream from them
		#bedtools slop -s -i $output/${sample}-$region.$nuc.top.tab -g $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes -b 3 > $output/$(basename $output/${sample}-$region.$nuc.top.tab .top.tab).flank.tab
	
		#Get nucleotide sequence of rNMPs and the 3 nucleotides up/downstream from them
		#bedtools getfasta -s -fi $fasta -bed $output/$(basename $output/${sample}-$region.$nuc.top.tab .top.tab).flank.tab | awk '/>/{$0 = ">" ++i substr($0, 2)} 1' - > $output/$(basename $output/${sample}-$region.$nuc.top.tab .top.tab).flank.txt

		#sites=$(grep -c "^>" $output/$(basename $output/${sample}-$region.$nuc.top.tab .top.tab).flank.txt)

		#meme $output/$(basename $output/${sample}-$region.$nuc.top.tab .top.tab).flank.txt -o $output/meme-$(basename $output/${sample}-$region.$nuc.top.tab .top.tab) -dna -minw 7 -nsites $sites -brief 1000000

		#done
	done
done

#rm $output/*.txt
