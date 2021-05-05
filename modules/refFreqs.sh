#!/usr/bin/env bash

#Author: Alli L. Gombolay
#Calculates nucleotide frequencies of reference genome

#############################################################################################################################

#Load config file
. "$1"

#############################################################################################################################

for file in $(dirname $fasta)/*.fa; do

	#region=$(basename $file .fa | cut -d "-" -f2)
		
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
	paste <(echo -e "A") <(echo "$A_BkgFreq" | xargs printf "%.*f\n" 5) > $(basename $file .fa).txt
	paste <(echo -e "C") <(echo "$C_BkgFreq" | xargs printf "%.*f\n" 5) >> $(basename $file .fa).txt
	paste <(echo -e "G") <(echo "$G_BkgFreq" | xargs printf "%.*f\n" 5) >> $(basename $file .fa).txt
	paste <(echo -e "U") <(echo "$T_BkgFreq" | xargs printf "%.*f\n" 5) >> $(basename $file .fa).txt

done
