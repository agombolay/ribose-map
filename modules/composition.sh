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
	
		other_new=$(echo "${other[*]}" | sed 's/ /|/g')
		
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito - | grep -Ewv ${other_new} -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -wv $mito $repository/results/$sample/coordinate$quality/$sample.bed | grep -Ewv ${other_new} - | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-nucleus.nucs.tab

		#Mito
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_mito.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_mito.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -w $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-mito.nucs.tab

		#Other
		for region in "${other[@]}"; do
			
			#Subset FASTA file based on region
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_$region.fa

			#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
			grep -w $region $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab
		
		done

	else
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito - )
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -wv $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-nucleus.nucs.tab

		#Mito
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_mito.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_mito.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -w $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-mito.nucs.tab
	fi
	
elif [[ ! $mito ]]; then
		
	if [[ $other ]]; then
		
		other_new=$(echo "${other[*]}" | sed 's/ /|/g')
		
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -Ewv ${other_new} -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa

		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		grep -Ewv ${other_new} $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-nucleus.nucs.tab

		#Other
		for region in "${other[@]}"; do
			
			#Subset FASTA file based on region
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_$region.fa

			#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
			grep -w $region $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab
		
		done

	else
		#Nucleus
		#Separate BED file by oraganelle and get nucleotide for each chromosomal coordinate
		bedtools getfasta -s -fi $fasta -bed $repository/results/$sample/coordinate$quality/$sample.bed | grep -v '>' > $output/${sample}-nucleus.nucs.tab
	fi

fi

######################################################################################################################################################
	
for file in $output/${sample}-*.nucs.tab; do

	temp=$(echo $file | awk -F '[-]' '{print $2 $3 $4}')
	region=$(basename $temp .nucs.tab)
	
	#Nucleotide Frequencies of Reference Genome
	
	#Calculate counts of each nucleotide
	A_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'A|a' - | wc -l)
	C_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'C|c' - | wc -l)
	G_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'G|g' - | wc -l)
	T_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'T|t' - | wc -l)
	
	#Calculate total number of nucleotides
	BkgTotal=$(($A_Bkg + $C_Bkg + $G_Bkg + $T_Bkg))
		
	#Calculate frequencies of each nucleotide
	A_BkgFreq=$(echo "($A_Bkg + $T_Bkg)/($BkgTotal*2)" | bc -l)
	C_BkgFreq=$(echo "($C_Bkg + $G_Bkg)/($BkgTotal*2)" | bc -l)
	G_BkgFreq=$(echo "($G_Bkg + $C_Bkg)/($BkgTotal*2)" | bc -l)
	T_BkgFreq=$(echo "($T_Bkg + $A_Bkg)/($BkgTotal*2)" | bc -l)
		
	#Save nucleotide frequencies to .txt file
	#echo $A_BkgFreq | xargs printf "%.*f\n" 5 > $output/A_Bkg.txt
	#echo $C_BkgFreq | xargs printf "%.*f\n" 5 > $output/C_Bkg.txt
	#echo $G_BkgFreq | xargs printf "%.*f\n" 5 > $output/G_Bkg.txt
	#echo $T_BkgFreq | xargs printf "%.*f\n" 5 > $output/T_Bkg.txt

	#Bkg=$(paste $output/{A,C,G,T}_Bkg.txt)

	#echo -e "A\tC\tG\tT" > "${fasta%.*}"-Freqs.$region.txt
	#paste <(echo -e "$Bkg") >> "${fasta%.*}"-Freqs.$region.txt
	
	paste <(echo -e "A") <(echo "$A_BkgFreq" | xargs printf "%.*f\n" 5) > "${fasta%.*}"-Freqs.$region.txt
	paste <(echo -e "C") <(echo "$C_BkgFreq" | xargs printf "%.*f\n" 5) >> "${fasta%.*}"-Freqs.$region.txt
	paste <(echo -e "G") <(echo "$G_BkgFreq" | xargs printf "%.*f\n" 5) >> "${fasta%.*}"-Freqs.$region.txt
	paste <(echo -e "U") <(echo "$T_BkgFreq" | xargs printf "%.*f\n" 5) >> "${fasta%.*}"-Freqs.$region.txt

######################################################################################################################################################
	
	#Nucleotide Frequencies of rNMPs
	if [ -s $output/${sample}-$region.nucs.tab ]; then
	
		#Calculate counts of each nucleotide
		A_Ribo=$(awk '$1 == "A" || $1 == "a"' $output/${sample}-$region.nucs.tab | wc -l)
		C_Ribo=$(awk '$1 == "C" || $1 == "c"' $output/${sample}-$region.nucs.tab | wc -l)
		G_Ribo=$(awk '$1 == "G" || $1 == "g"' $output/${sample}-$region.nucs.tab | wc -l)
		U_Ribo=$(awk '$1 == "T" || $1 == "t"' $output/${sample}-$region.nucs.tab | wc -l)
	
		#Calculate total number of nucleotides
		RiboTotal=$(($A_Ribo + $C_Ribo + $G_Ribo + $U_Ribo))
	
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
	
		#Save nucleotide frequencies to .txt file
		paste <(echo -e "rA") <(echo "$A_RiboFreq2" | xargs printf "%.*f\n" 5) > $output/${sample}-$region.counts.txt
		paste <(echo -e "rC") <(echo "$C_RiboFreq2" | xargs printf "%.*f\n" 5) >> $output/${sample}-$region.counts.txt
		paste <(echo -e "rG") <(echo "$G_RiboFreq2" | xargs printf "%.*f\n" 5) >> $output/${sample}-$region.counts.txt
		paste <(echo -e "rU") <(echo "$U_RiboFreq2" | xargs printf "%.*f\n" 5) >> $output/${sample}-$region.counts.txt
	fi
done

######################################################################################################################################################

#Remove temporary files
#rm $output/*.nucs.tab

#Print status
#echo "Status: Composition Module for $sample is complete"
