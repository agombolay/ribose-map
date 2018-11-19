#!/usr/bin/env bash

#Alli Gombolay, 11/2018
#Calculate top 1% of sites: https://www.dummies.com/education/math/statistics/how-to-calculate-percentiles-in-statistics/

######################################################################################################################################################

directory=$1
output=$directory/percentiles; rm -rf $output; mkdir $output

######################################################################################################################################################

for file in $(ls $directory/counts/*.tab); do

	sample=$(basename $file | cut -d. -f1)

	#Calculate index
	index=$(echo "$(wc -l < $file)*0.99" | bc)

	#Test if index is integer or floating
	integer=$(echo $index | grep -E '[0-9]+\.[0]{2}' -)
	floating=$(echo $index | grep -E '[0-9]+\.[^00]' -)

	#Sort from low to highest rNMP counts
	sort -k5,5n $file > $output/$sample.sorted.tsv

######################################################################################################################################################
	if [[ $integer ]]; then
		
		#Calculate line numbers
		line1=$(echo $index | awk '{print int($1)}')
		line2=$(echo "$line1+1" | bc)

		#Get counts of rNMPs for lines one and two
		value1=$(head -${line1} $output/$sample.sorted.tsv | tail -1 | awk '{print $5}')
		value2=$(head -${line2} $output/$sample.sorted.tsv | tail -1 | awk '{print $5}')

		#Calculate 99th percentile (average of counts)
		percentile=$(echo "scale=2; ($value1 + $value2)/2" | bc -l)

	elif [[ $floating ]]; then
		
		#Round up value of index to get line number
		line=$(echo $index | awk '{print int($1 + 1)}')

		#Calculate 99th percentile (count at the line)
		percentile=$(head -${line} $output/$sample.sorted.tsv | tail -1 | awk '{print $5}')
		
	fi
######################################################################################################################################################

	#Save top 1% of data to new file
	mawk -v "OFS=\t" -v "x=$percentile" '{if ($5 >= x) print $0}' $output/$sample.sorted.tsv | sort > $output/$sample.top1%.tsv

done

#Remove temporary files
rm -f $output/*.sorted.tsv
