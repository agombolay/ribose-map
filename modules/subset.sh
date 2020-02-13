#!/usr/bin/env bash

#Load config file
. "$1"

#Specify directory
output=$repository/results/$sample/coordinate$quality

if [[ $mito ]]; then
	
	#Mito
	grep -w $mito $output/$sample.bed > $output/${sample}-$mito.coords.bed
	cut -f1,2,3,6 $output/${sample}-$mito.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-mito.counts.tab

	chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)

	samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$mito.fa
	samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$mito.fa

	if [[ $nucleus ]]; then
	
		if [[ $other ]]; then

			other_new=$(echo $other | sed 's/ /|/g')
		
			#Nucleus
			grep -wv $mito $output/$sample.bed | grep -Ewv $other_new - > $output/${sample}-nucleus.coords.bed
			cut -f1,2,3,6 $output/${sample}-nucleus.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-nucleus.counts.tab

			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito - | grep -Ewv $other_new -)

			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

			#Other
			for region in $other; do
				
				grep -w $region $output/$sample.bed > $output/${sample}-$region.coords.bed
				cut -f1,2,3,6 $output/${sample}-$region.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-$region.counts.tab
			
				chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)

				samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
				samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

			done

		else

			#Nucleus
			grep -wv $mito $output/$sample.bed > $output/${sample}-nucleus.coords.bed
			cut -f1,2,3,6 $output/${sample}-nucleus.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-nucleus.counts.tab
		
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito -)

			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

		fi

	else

		#Other
		for region in $other; do

			grep -w $region $output/$sample.bed > $output/${sample}-$region.coords.bed
			cut -f1,2,3,6 $output/${sample}-$region.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-$region.counts.tab
		
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)

			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

		done

	fi
	
else

	if [[ $nucleus ]]; then
		
		if [[ $other ]]; then
		
			other_new=$(echo $other | sed 's/ /|/g')
		
			#Nucleus
			grep -Ewv $other_new $output/$sample.bed > $output/${sample}-nucleus.coords.bed
			cut -f1,2,3,6 $output/${sample}-nucleus.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-nucleus.counts.tab
		
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -Ewv $other_new -)

			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

			#Other
			for region in $other; do

				grep -w $region $output/$sample.bed > $output/${sample}-$region.coords.bed
				cut -f1,2,3,6 $output/${sample}-$region.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-$region.counts.tab
				
				chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)

				samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
				samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

			done

		else

			#Nucleus
			cp $output/$sample.bed $output/$sample-nucleus.coords.bed
			cut -f1,2,3,6 $output/${sample}-nucleus.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-nucleus.counts.tab
		
			cp $fasta $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

		fi

	else

		#Other
		for region in $other; do

			grep -w $region $output/$sample.bed > $output/${sample}-$region.coords.bed
			cut -f1,2,3,6 $output/${sample}-$region.coords.bed | uniq -c - | awk -v "OFS=\t" '{print $2, $3, $4, $5, $1}' - | sort -k5,5n - > $output/${sample}-$region.counts.tab
		
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)

			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

		done

	fi
fi
