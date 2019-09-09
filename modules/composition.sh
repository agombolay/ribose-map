#!/usr/bin/env bash

#Alli Gombolay, 11/2018
#Calculate counts of rAMP, rCMP, rGMP, and rUMP

######################################################################################################################################################
#Load config file
. "$1"

output=$repository/results/$sample/composition$quality
rm -r $output; mkdir -p $output

######################################################################################################################################################

for region in "nucleus" "$mito" "$other"; do

	#Separate BED file by oraganelle
	if [[ $region == "nucleus" ]]; then
		#Get nucleotide for each chromosomal coordinate			
		grep -wv $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab

	elif [[ $region == "$mito" ]]; then
		#Get nucleotide for each chromosomal coordinate
		grep -w $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab

	elif [[ $region == "$other" ]]; then
		#Get nucleotide for each chromosomal coordinate
		grep -w $other $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab

	fi

	A_Ribo=$(awk '$1 == "A" || $1 == "a"' $output/${sample}-$region.nucs.tab | wc -l)
	C_Ribo=$(awk '$1 == "C" || $1 == "c"' $output/${sample}-$region.nucs.tab | wc -l)
	G_Ribo=$(awk '$1 == "G" || $1 == "g"' $output/${sample}-$region.nucs.tab | wc -l)
	U_Ribo=$(awk '$1 == "T" || $1 == "t"' $output/${sample}-$region.nucs.tab | wc -l)
	
	RiboTotal=$(($A_Ribo + $C_Ribo + $G_Ribo + $U_Ribo))

	A_RiboFreq=$(echo "($A_Ribo)" | bc -l)
	C_RiboFreq=$(echo "($C_Ribo)" | bc -l)
	G_RiboFreq=$(echo "($G_Ribo)" | bc -l)
	U_RiboFreq=$(echo "($U_Ribo)" | bc -l)
	
	paste <(echo -e "rA") <(echo "$A_RiboFreq") >> $output/${sample}-$region.counts.tab
	paste <(echo -e "rC") <(echo "$C_RiboFreq") >> $output/${sample}-$region.counts.tab
	paste <(echo -e "rG") <(echo "$G_RiboFreq") >> $output/${sample}-$region.counts.tab
	paste <(echo -e "rU") <(echo "$U_RiboFreq") >> $output/${sample}-$region.counts.tab

done

######################################################################################################################################################

#Remove temporary files
rm $output/*.nucs.tab

#Print status
echo "Status: Composition Module for $sample is complete"
