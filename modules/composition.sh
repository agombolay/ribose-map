#!/usr/bin/env bash

#Author: Alli L. Gombolay
#Calculates normalized percentage of each rNMP

######################################################################################################################################################

#Load config file
. "$1"

output=$repository/results/$sample/composition$quality
rm -r $output; mkdir -p $output

######################################################################################################################################################

for region in $units "chromosomes"; do

	#Get nucleotide for each genomic coordinate
	bedtools getfasta -s -fi $fasta -bed $repository/results/$sample/coordinate$quality/${sample}-$region.coords.bed | grep -v '>' > $output/${sample}-$region.nucs.tab
	
	#Calculate counts of each nucleotide
	A_Ribo=$(awk '$1 == "A" || $1 == "a"' $output/${sample}-$region.nucs.tab | wc -l)
	C_Ribo=$(awk '$1 == "C" || $1 == "c"' $output/${sample}-$region.nucs.tab | wc -l)
	G_Ribo=$(awk '$1 == "G" || $1 == "g"' $output/${sample}-$region.nucs.tab | wc -l)
	U_Ribo=$(awk '$1 == "T" || $1 == "t"' $output/${sample}-$region.nucs.tab | wc -l)
	
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

done

######################################################################################################################################################

#Remove temporary files
rm $output/*.nucs.tab

#Print status
echo "Status: Composition Module for $sample is complete"
