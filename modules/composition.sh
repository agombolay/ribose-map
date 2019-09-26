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
		grep -wv $mito $repository/results/$sample/coordinate$quality/$sample.bed | grep -wv $other - | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab

	elif [[ $region == "$mito" ]]; then
		#Get nucleotide for each chromosomal coordinate
		grep -w $mito $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab

	elif [[ $region == "$other" ]]; then
		#Get nucleotide for each chromosomal coordinate
		grep -w $other $repository/results/$sample/coordinate$quality/$sample.bed | bedtools getfasta -s -fi $fasta -bed - | grep -v '>' > $output/${sample}-$region.nucs.tab

	fi

	if [[ -s $output/${sample}-$region.nucs.tab ]]; then
	
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
	fi
done

######################################################################################################################################################

#Create .fai file for reference
samtools faidx $fasta

#Create .bed file for reference
cut -f 1,2 $fasta.fai > $(dirname $fasta)/$(basename $fasta .fa).bed

######################################################################################################################################################

#Subset FASTA file based on region
for region in "nucleus" "$mito" "$other; do

	if [[ $region == "nucleus" ]]; then
			
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito - | grep -wv $other -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_nucleus.fa
		fi
		
	elif [[ $region == "mitochondria" ]]; then
		
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_mitochondria.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_mitochondria.fa
		fi

	elif [[ $region == "other" ]]; then
		
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $other -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)_mitochondria.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)_mitochondria.fa
		fi
	fi
	
######################################################################################################################################################

	#Calculate counts of each nucleotide
	A_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'A|a' - | wc -l)
	C_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'C|c' - | wc -l)
	G_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'G|g' - | wc -l)
	T_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)_$region.fa | grep -Eo 'T|t' - | wc -l)
	
	#Calculate total number of nucleotides
	BkgTotal=$(($A_Bkg + $C_Bkg + $G_Bkg + $T_Bkg))
		
	#Calculate frequency of each nucleotide
	A_BkgFreq=$(echo "($A_Bkg + $T_Bkg)/($BkgTotal*2)" | bc -l)
	C_BkgFreq=$(echo "($C_Bkg + $G_Bkg)/($BkgTotal*2)" | bc -l)
	G_BkgFreq=$(echo "($G_Bkg + $C_Bkg)/($BkgTotal*2)" | bc -l)
	T_BkgFreq=$(echo "($T_Bkg + $A_Bkg)/($BkgTotal*2)" | bc -l)
		
	#Save background frequencies of dNMPs to TXT files
	echo $A_BkgFreq | xargs printf "%.*f\n" 5 > $output/A_Bkg.txt
	echo $C_BkgFreq | xargs printf "%.*f\n" 5 > $output/C_Bkg.txt
	echo $G_BkgFreq | xargs printf "%.*f\n" 5 > $output/G_Bkg.txt
	echo $T_BkgFreq | xargs printf "%.*f\n" 5 > $output/T_Bkg.txt

	#Combine dNMP frequencies into one file
	Bkg=$(paste $output/{A,C,G,T}_Bkg.txt)

done

######################################################################################################################################################

#Remove temporary files
rm $output/*.nucs.tab

#Print status
echo "Status: Composition Module for $sample is complete"
