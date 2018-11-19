#!/usr/bin/env bash

#Alli Gombolay, 11/2018
#Create lists of hotspots

######################################################################################################################################################

while getopts ":h" opt;
    do case ${opt} in
        h ) echo "Usage: ./hotspots.sh /path/to/BED /path/to/FASTA" ;;
        \? ) echo "For usage statement: ./hotspots.sh -h" ;;
    esac
done

directory=$1
reference=$2
output=$directory/hotspots; rm -rf $output; mkdir $output

######################################################################################################################################################

/users/alligombolay/desktop/hotspots/counts.sh $directory $reference

for option in "all" "percentile"; do

######################################################################################################################################################
	if [[ $option == "all" ]]; then
			
		#CREATE FILE 1
		#Combine files into one
		cat $directory/counts/*.tab >> $output/combined.unsorted.$option.tsv

		#Print total number of libraries and number of libraries that have a particular coordinate
		sort -k1,1 -k2,2n -k4,4 -k5 $output/combined.unsorted.$option.tsv | cut -f1,2,3,4,5,6 | uniq -c | awk -v "OFS=\t" \
		-v "x=$(ls $directory/counts/*.tab | wc -l)" '{print $2 "_" $3 "_" $4 "_" $5, $6, $7, $1, x}' > $output/combined.sorted.$option.tsv

		#CREATE FILE 2
		#Combine coordinates into one column using the underscore character and add name of library
		for file in $(ls $directory/counts/*.tab); do

			sample=$(basename "$file" | cut -d. -f1)
			cut -f1,2,3,4 $file | awk -v "OFS=\t" -v "y=$sample" '{print $1 "_" $2 "_" $3 "_" $4, y}' | sort > $output/$sample.temp1.$option.tsv

		done

		#Join first two temporary files together based on their chromosomal coordinates
		file1=$(ls $output/*.temp1.$option.tsv | head -1)
		file2=$(ls $output/*.temp1.$option.tsv | head -2 | tail -1)

		command="join -t $'\t' $file1 $file2 -a1 -a2"
		eval $command

		#Join previous output with 3rd, 4th files etc. until all files have been joined
		for file in $(ls $output/*.temp1.$option.tsv | tail -n+3); do
			command+="| join -t $'\t' $file - -a1 -a2"
		done

		command+="> $output/combined.joined.$option.tsv"
		eval $command

		#Join two files together into one
		join -t $'\t' <(sort $output/combined.sorted.$option.tsv) <(sort $output/combined.joined.$option.tsv) | sed 's/_/\t/g' | sort -nrk 5,5 \
		> $output/hotspots.$option.tsv

######################################################################################################################################################
	elif [[ $option == "percentile" ]]; then
			
		#CREATE FILE 1
		#Combine files into one
		/users/alligombolay/desktop/hotspots/percentile.sh $directory
		cat $directory/percentiles/*top1%.tsv >> $output/combined.unsorted.$option.tsv

		#Print total number of libraries and number of libraries that have a particular coordinate
		sort -k1,1 -k2,2n -k4,4 -k5 $output/combined.unsorted.$option.tsv | cut -f1,2,3,4,5,6 | uniq -c | awk -v "OFS=\t" \
		-v "x=$(ls $directory/counts/*.tab | wc -l)" '{print $2 "_" $3 "_" $4 "_" $5, $6, $7, $1, x}' > $output/combined.sorted.$option.tsv

		#CREATE FILE 2
		#Combine coordinates into one column using the underscore character and add name of library
		for file in $(ls $directory/percentiles/*top1%.tsv); do

			sample=$(basename "$file" | cut -d. -f1)
			cut -f1,2,3,4 $file | awk -v "OFS=\t" -v "y=$sample" '{print $1 "_" $2 "_" $3 "_" $4, y}' | sort > $output/$sample.temp1.$option.tsv

		done

		#Join first two temporary files together based on their chromosomal coordinates
		file1=$(ls $output/*.temp1.$option.tsv | head -1)
		file2=$(ls $output/*.temp1.$option.tsv | head -2 | tail -1)

		command="join -t $'\t' $file1 $file2 -a1 -a2"
		eval $command

		#Join previous output with 3rd, 4th files etc. until all files have been joined
		for file in $(ls $output/*.temp1.$option.tsv | tail -n+3); do
			command+="| join -t $'\t' $file - -a1 -a2"
		done

		command+="> $output/combined.joined.$option.tsv"
		eval $command

		#Join two files together into one
		join -t $'\t' <(sort $output/combined.sorted.$option.tsv) <(sort $output/combined.joined.$option.tsv) | sed 's/_/\t/g' | sort -nrk 5,5 \
		> $output/hotspots.$option.tsv

	fi
######################################################################################################################################################
done

#Remove temporary files
rm -f $output/combined.*.tsv $output/*.temp1.*.tsv

#Print status
echo "Status: Hotspot Module is complete"
