#!/usr/bin/env bash

#Alli Gombolay, 11/2018
#Calculate counts of rAMP, rCMP, rGMP, and rUMP

######################################################################################################################################################
#Load config file
. "$1"

output=$repository/results/$sample/composition$quality
rm -r $output; mkdir -p $output

######################################################################################################################################################

#Create .fai file for reference
samtools faidx $fasta

#Create .bed file for reference
cut -f 1,2 $fasta.fai > $(dirname $fasta)/$(basename $fasta .fa).bed

######################################################################################################################################################

if [[ $mito ]]; then
	
	if [[ $other ]]; then
	
		other_new=$(echo $other | sed 's/ /|/g')
		
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito - | grep -Ewv $other_new -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -wv $mito $repository/results/$sample/coordinate$quality/$sample.bed | grep -Ewv $other_new - | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-nucleus.nucs.tab

		#Mito
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-mito.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-mito.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -w $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-mito.nucs.tab

		#Other
		for region in $other; do
			
			#Subset FASTA file based on region
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

			#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
			grep -w $region $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab
		
		done

	else
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito - )
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -wv $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-nucleus.nucs.tab

		#Mito
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-mito.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-mito.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -w $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-mito.nucs.tab
	fi
	
elif [[ ! $mito ]]; then
		
	if [[ $other ]]; then
		
		other_new=$(echo $other | sed 's/ /|/g')
		
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -Ewv $other_new -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -Ewv $other_new $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-nucleus.nucs.tab

		#Other
		for region in $other; do
			
			#Subset FASTA file based on region
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

			#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
			grep -w $region $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab
		
		done

	else
		#Nucleus
		#Subset FASTA file based on region
		cat $fasta > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
		
		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		bedtools getfasta -s -fi $fasta -bed $repository/results/$sample/coordinate$quality/$sample.bed | grep -v '>' > $output/${sample}-nucleus.nucs.tab
	fi

fi

######################################################################################################################################################
	
for file in $(dirname $fasta)/$(basename $fasta .fa)-*.fa; do

	temp=$(echo $file | awk -F '[-]' '{print $2 $3 $4}')
	region=$(basename $temp .fa)
	
	#Nucleotide Frequencies of Reference Genome
	
	#Calculate counts of each nucleotide
	A_Bkg=$(grep -v '>' $file | grep -Eo 'A|a' - | wc -l)
	C_Bkg=$(grep -v '>' $file | grep -Eo 'C|c' - | wc -l)
	G_Bkg=$(grep -v '>' $file | grep -Eo 'G|g' - | wc -l)
	T_Bkg=$(grep -v '>' $file | grep -Eo 'T|t' - | wc -l)
	
	#Calculate total number of nucleotides
	BkgTotal=$(($A_Bkg + $C_Bkg + $G_Bkg + $T_Bkg))
		
	#Calculate frequencies of each nucleotide
	A_BkgFreq=$(echo "($A_Bkg + $T_Bkg)/($BkgTotal*2)" | bc -l)
	C_BkgFreq=$(echo "($C_Bkg + $G_Bkg)/($BkgTotal*2)" | bc -l)
	G_BkgFreq=$(echo "($G_Bkg + $C_Bkg)/($BkgTotal*2)" | bc -l)
	T_BkgFreq=$(echo "($T_Bkg + $A_Bkg)/($BkgTotal*2)" | bc -l)
		
	#Save nucleotide frequencies to .txt file
	paste <(echo -e "A") <(echo "$A_BkgFreq" | xargs printf "%.*f\n" 5) > $(dirname $fasta)/$(basename $fasta .fa)-$region.txt
	paste <(echo -e "C") <(echo "$C_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.txt
	paste <(echo -e "G") <(echo "$G_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.txt
	paste <(echo -e "U") <(echo "$T_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.txt

done

######################################################################################################################################################

for file in $output/${sample}-*.nucs.tab; do

	#Nucleotide Frequencies of rNMPs
	if [ -s $output/${sample}-*.nucs.tab ]; then
	
		temp=$(echo $file | awk -F '[-]' '{print $2 $3 $4}')
		region=$(basename $temp .nucs.tab)
	
		#Calculate counts of each nucleotide
		A_Ribo=$(awk '$1 == "A" || $1 == "a"' $file | wc -l)
		C_Ribo=$(awk '$1 == "C" || $1 == "c"' $file | wc -l)
		G_Ribo=$(awk '$1 == "G" || $1 == "g"' $file | wc -l)
		U_Ribo=$(awk '$1 == "T" || $1 == "t"' $file | wc -l)
	
		#Calculate total number of nucleotides
		RiboTotal=$(($A_Ribo + $C_Ribo + $G_Ribo + $U_Ribo))
	
		A_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -1)
		C_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -2 | tail -1)
		G_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -3 | tail -1)
		T_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -4 | tail -1)
	
		#Calculate normalized frequency of each nucleotide, step 1
		A_RiboFreq1=$(echo "($A_Ribo/$RiboTotal)/$A_BkgFreq" | bc -l)
		C_RiboFreq1=$(echo "($C_Ribo/$RiboTotal)/$C_BkgFreq" | bc -l)
		G_RiboFreq1=$(echo "($G_Ribo/$RiboTotal)/$G_BkgFreq" | bc -l)
		U_RiboFreq1=$(echo "($U_Ribo/$RiboTotal)/$T_BkgFreq" | bc -l)
	
		#Calculate normalized frequency of each nucleotide, step 2
		A_RiboFreq2=$(echo "$A_RiboFreq1/($A_RiboFreq1 + $C_RiboFreq1 + $G_RiboFreq1 + $U_RiboFreq1)*100" | bc -l)
		C_RiboFreq2=$(echo "$C_RiboFreq1/($A_RiboFreq1 + $C_RiboFreq1 + $G_RiboFreq1 + $U_RiboFreq1)*100" | bc -l)
		G_RiboFreq2=$(echo "$G_RiboFreq1/($A_RiboFreq1 + $C_RiboFreq1 + $G_RiboFreq1 + $U_RiboFreq1)*100" | bc -l)
		U_RiboFreq2=$(echo "$U_RiboFreq1/($A_RiboFreq1 + $C_RiboFreq1 + $G_RiboFreq1 + $U_RiboFreq1)*100" | bc -l)
	
		#Save nucleotide counts to .txt file
		paste <(echo -e "rA") <(echo "$A_Ribo") > $output/${sample}-$region.counts.txt
		paste <(echo -e "rC") <(echo "$C_Ribo") >> $output/${sample}-$region.counts.txt
		paste <(echo -e "rG") <(echo "$G_Ribo") >> $output/${sample}-$region.counts.txt
		paste <(echo -e "rU") <(echo "$U_Ribo") >> $output/${sample}-$region.counts.txt
		
		#Save nucleotide frequencies to .txt file
		paste <(echo -e "rA") <(echo "$A_RiboFreq2" | xargs printf "%.*f\n" 5) > $output/${sample}-$region.frequencies.txt
		paste <(echo -e "rC") <(echo "$C_RiboFreq2" | xargs printf "%.*f\n" 5) >> $output/${sample}-$region.frequencies.txt
		paste <(echo -e "rG") <(echo "$G_RiboFreq2" | xargs printf "%.*f\n" 5) >> $output/${sample}-$region.frequencies.txt
		paste <(echo -e "rU") <(echo "$U_RiboFreq2" | xargs printf "%.*f\n" 5) >> $output/${sample}-$region.frequencies.txt
	
	fi
done

######################################################################################################################################################

#Remove temporary files
rm $output/*.nucs.tab

#Print status
echo "Status: Composition Module for $sample is complete"
