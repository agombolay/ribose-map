#!/usr/bin/env bash

#Alli Gombolay, 02/2020
#Subsets data by genomic region

######################################################################################################################################################
#Load config file
. "$1"

#Specify directory
output=$repository/results/$sample/coordinate$quality

#Create index file
samtools faidx $fasta

######################################################################################################################################################

if [[ $other ]]; then

	other_new=$(echo $other | sed 's/ /|/g')
		
	#Chromosomes
	grep -Ewv $other_new $output/$sample.bed > $output/${sample}-chromosomes.coords.bed
	cut -f1,2,3,6 $output/${sample}-chromosomes.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-chromosomes.counts.tab

	#Create FASTA and FAI files for Chromosomes
	chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -Ewv $other_new -)

	samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa
	samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa
	
	#Create BED file for Chromosomes
	cut -f 1,2 $fasta.fai > $(dirname $fasta)/$(basename $fasta .fa).bed

	#Other
	for region in $other; do
				
		grep -w $region $output/$sample.bed > $output/${sample}-$region.coords.bed
		cut -f1,2,3,6 $output/${sample}-$region.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-$region.counts.tab
		
		#Create FASTA and FAI files for Other
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)

		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
		
		#Create BED file for Other
		cut -f 1,2 $fasta.fai > $(dirname $fasta)/$(basename $fasta .fa).bed

	done
else
	
	#Chromosomes
	cut -f1,2,3,6 $output/${sample}.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-chromosomes.counts.tab

fi
