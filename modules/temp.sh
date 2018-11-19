#!/usr/bin/env bash

#Alli Gombolay, 11/2018
#Create file of ribo counts with flanking nucleotides

directory=$1
reference=$2

output=$directory/counts; rm -rf $output; mkdir $output

######################################################################################################################################################

#Create FASTA index file and BED file for reference
samtools faidx $reference && cut -f 1,2 $reference.fai > $output/reference.bed

######################################################################################################################################################

for file in $(ls $directory/*.bed); do 

	sample=$(basename "$file" | cut -d. -f1)

	#Separate BED file by oraganelle
	if [ $region == "nucleus"]; then
		#Get nucleotide for each chromosomal coordinate			
		grep -wvE 'chrM' $file > $output/$sample.$region.bed > $output/$sample.$region.bed

	elif [ $region == "mitochondria"]; then
		#Get nucleotide for each chromosomal coordinate
		grep -wE 'chrM' $file > $output/$sample.$region.bed > $output/$sample.$region.bed
	fi


######################################################################################################################################################

	for region in "nucleus" "mitochondria"; do
		
		#Get nucleotide at each ribonucleotide
		bedtools getfasta -s -fi $reference -bed $output/$sample.$region.bed | grep -v '>' > $output/$sample.$region.ribos.txt

		#Get nucleotides 2 bp upstream and downstream
		bedtools slop -s -i $output/$sample.$region.bed -g $output/reference.bed -b 2 | bedtools getfasta -s -fi $reference -bed - | \
		grep -v '>' > $output/$sample.$region.flank.txt
	
		#Get counts of each ribonucleotide coordinate
		paste $file $output/$sample.$region.ribos.txt $output/$sample.$region.flank.txt | cut -f1,2,3,6,7,8 - | uniq -c - | \
		mawk -v "OFS=\t" '{print $2, $3, $4, $5, $6, $7, $1}' > $output/$sample.$region.counts.tab

	done

done

######################################################################################################################################################

#Remove temporary files
rm $output/*.{txt, bed}
