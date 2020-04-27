#!/usr/bin/env bash

#Alli Gombolay, 02/2020
#Subsets data by genomic region

######################################################################################################################################################

#Load config file
. "$1"

######################################################################################################################################################

if [[ $other ]]; then

	other_new=$(echo $other | sed 's/ /|/g')
		
	#Chromosomes
	grep -Ewv $other_new $repository/results/$sample/coordinate$quality/$sample.bed > $repository/results/$sample/coordinate$quality/${sample}-chromosomes.coords.bed
	cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-chromosomes.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-chromosomes.counts.tab

	#Create FASTA and FAI files for Chromosomes
	chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes | grep -Ewv $other_new -)

	samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa
	samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa

	#Other
	for region in $other; do
				
		grep -w $region $repository/results/$sample/coordinate$quality/$sample.bed > $repository/results/$sample/coordinate$quality/${sample}-$region.coords.bed
		cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-$region.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-$region.counts.tab
		
		#Create FASTA and FAI files for Other
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes | grep -w $region -)

		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

	done
else
	
	#Chromosomes
	cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-chromosomes.counts.tab

fi
