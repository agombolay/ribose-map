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
	grep -Ewv $other_new $repository/results/$sample/coordinate$quality/$sample.bed > $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed
	cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-chromosomes.Combined.tab

	#Create FASTA and FAI files for Chromosomes
	chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes | grep -Ewv $other_new -)

	samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa
	samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa

	for nuc in "A" "C" "G" "U"; do

		#Save frequencies of rNMPs to TXT files
		if [[ $nuc == "A" ]]; then
			#Create BED file for only A ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "A" || $2 == "a"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab
			
		elif [[ $nuc == "C" ]]; then
			#Create BED file for only C ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "C" || $2 == "c"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab
			
		elif [[ $nuc == "G" ]]; then	
			#Create BED file for only G ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "G" || $2 == "g"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab
			
		elif [[ $nuc == "U" ]]; then	
			#Create BED file for only U ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "T" || $2 == "t"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab
		fi
	done
	
	#Other
	for region in $other; do
				
		grep -w $region $repository/results/$sample/coordinate$quality/$sample.bed > $repository/results/$sample/coordinate$quality/${sample}-$region.bed
		cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-$region.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-$region.Combined.tab
		
		#Create FASTA and FAI files for Other
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes | grep -w $region -)

		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
		
		for nuc in "A" "C" "G" "U"; do

			#Save frequencies of rNMPs to TXT files
			if [[ $nuc == "A" ]]; then
				#Create BED file for only A ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.bed | awk '$2 == "A" || $2 == "a"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed
				cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.tab

			elif [[ $nuc == "C" ]]; then
				#Create BED file for only C ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.bed | awk '$2 == "C" || $2 == "c"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed
				cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.tab

			elif [[ $nuc == "G" ]]; then	
				#Create BED file for only G ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.bed | awk '$2 == "G" || $2 == "g"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed
				cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.tab

			elif [[ $nuc == "U" ]]; then	
				#Create BED file for only U ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.bed | awk '$2 == "T" || $2 == "t"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed
				cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.$region.$nuc.tab
			fi
		done

	done
else
	
	#Chromosomes
	cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-chromosomes.Combined.tab

	for nuc in "A" "C" "G" "U"; do

		#Save frequencies of rNMPs to TXT files
		if [[ $nuc == "A" ]]; then
			#Create BED file for only A ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "A" || $2 == "a"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab

		elif [[ $nuc == "C" ]]; then
			#Create BED file for only C ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "C" || $2 == "c"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab

		elif [[ $nuc == "G" ]]; then	
			#Create BED file for only G ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "G" || $2 == "g"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab

		elif [[ $nuc == "U" ]]; then	
			#Create BED file for only U ribonucleotide
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-chromosomes.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-chromosomes.bed | awk '$2 == "T" || $2 == "t"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed
			cut -f1,2,3,6 $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, ".", ".", $5, $1}' - | sort -k7,7n - > $repository/results/$sample/coordinate$quality/${sample}-Ribo.chromosomes.$nuc.tab
		fi
	done
fi
