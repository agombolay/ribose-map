#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates rNMP frequencies and flanking dNMP frequencies (+/- 100 bp)

#Usage statement
function usage () {
	echo "Usage: Frequencies.sh [-s] 'Sample(s)' [-r] 'Reference' [-d] 'Directory' [-h]
	-s Sample name(s) (e.g., FS1, FS2, FS3 etc.)
	-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
	-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
        s ) sample=($OPTARG) ;;
	#Allow only one input argument
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Calculate frequencies
for sample in ${sample[@]}; do
	for subset in "all"; do
	
#############################################################################################################################
	#Input files
	reads=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-ReadInformation.$subset.txt
	BED=$directory/Ribose-Map/Reference/$reference.bed; FASTA=$directory/Ribose-Map/Reference/$reference.fa
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-Coordinates.$subset.bed

	#Output directory and file
	output=$directory/Ribose-Map/Results/$reference/$sample/Frequencies
	dataset=$output/$sample-Frequencies.$reference.$subset.txt
	
	#Create directory and remove old file
	mkdir -p $output; rm -f $output/$dataset
	
#############################################################################################################################
	#STEP 1: Calculate frequencies of reference genome
	
	#Subset FASTA file based on region
	if [ $subset == "all" ]; then
		cat $FASTA > temp.fa; samtools faidx temp.fa
	elif [ $subset == "mito" ]; then
		chromosomes=$(awk '{print $1}' $BED | grep -E '(chrM|MT)')
		samtools faidx $FASTA $chromosomes > temp.fa; samtools faidx temp.fa
	elif [ $subset == "nucleus" ]; then
		chromosomes=$(awk '{print $1}' $BED | grep -vE '(chrM|MT)')
		samtools faidx $FASTA $chromosomes > temp.fa; samtools faidx temp.fa
	fi

	#Calculate counts of each nucleotide
	A_BkgCount=$(grep -v '>' temp.fa | grep -o 'A' - | wc -l); C_BkgCount=$(grep -v '>' temp.fa | grep -o 'C' - | wc -l)
	G_BkgCount=$(grep -v '>' temp.fa | grep -o 'G' - | wc -l); T_BkgCount=$(grep -v '>' temp.fa | grep -o 'T' - | wc -l)
	
	#Calculate total number of nucleotides
	total_Bkg=$(($A_BkgCount+$C_BkgCount+$G_BkgCount+$T_BkgCount))

	#Calculate frequency of each nucleotide
	A_BkgFreq=$(echo "scale=12; $A_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')
	C_BkgFreq=$(echo "scale=12; $C_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')
	G_BkgFreq=$(echo "scale=12; $G_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')
	T_BkgFreq=$(echo "scale=12; $T_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')
	
#############################################################################################################################
	#STEP 2: Calculate frequencies of rNMPs in libraries

	if [ $subset == "all" ]; then
		#Select all reads located and extract rNMP bases
		cat $reads | awk '{print substr($0,length($0))}' - > RiboBases.txt
	elif [ $subset == "mito" ]; then
		#Select only reads located in mito DNA and extract rNMP bases
		grep -E '(chrM|MT)' $reads | awk '{print substr($0,length($0))}' - > RiboBases.txt
	elif [ $subset == "nucleus" ]; then
		#Select only reads located in nuclear DNA and extract rNMP bases
		grep -v -E '(chrM|MT)' $reads | awk '{print substr($0,length($0))}' - > RiboBases.txt
	fi
	
	#Calculate counts of rNMPs
	A_RiboCount=$(awk '$1 == "A"' RiboBases.txt | wc -l); C_RiboCount=$(awk '$1 == "C"' RiboBases.txt | wc -l)
	G_RiboCount=$(awk '$1 == "G"' RiboBases.txt | wc -l); U_RiboCount=$(awk '$1 == "T"' RiboBases.txt | wc -l)
	
	#Calculate total number of rNMPs
	RiboCount=$(($A_RiboCount+$C_RiboCount+$G_RiboCount+$U_RiboCount))

	#Calculate normalized frequency of each rNMP
	A_RiboFreq=$(echo "scale=12; ($A_RiboCount/$RiboCount)/$A_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	C_RiboFreq=$(echo "scale=12; ($C_RiboCount/$RiboCount)/$C_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	G_RiboFreq=$(echo "scale=12; ($G_RiboCount/$RiboCount)/$G_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	U_RiboFreq=$(echo "scale=12; ($U_RiboCount/$RiboCount)/$T_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')

	#Save normalized frequencies of rNMPs together
	RiboFreq=$(echo -e "$A_RiboFreq\t$C_RiboFreq\t$G_RiboFreq\t$U_RiboFreq")
	
#############################################################################################################################
	#STEP 3: Obtain coordinates/sequences of dNMPs +/- 100 bp from rNMPs

	#Obtain coordinates of sequences upstream/downstream from rNMPs
	bedtools flank -i $coordinates -s -g $BED -l 100 -r 0 | awk '$2 != "0" && $3 != "0"' - > Upstream.bed
	bedtools flank -i $coordinates -s -g $BED -l 0 -r 100 | awk '$2 != "0" && $3 != "0"' - > Downstream.bed
	
	#Obtain sequences of the sequences upstream/downstream from rNMPs
	bedtools getfasta -s -fi temp.fa -bed Upstream.bed -fo Upstream.fa
	bedtools getfasta -s -fi temp.fa -bed Downstream.bed -fo Downstream.fa
	
#############################################################################################################################
	#STEP 4: Insert tabs between sequences of dNMPs +/- 100 bp from rNMPs

	#Extract sequences from FASTA files (reverse order of upstream) and insert tabs between each nucleotide
	grep -v '>' Upstream.fa | rev > Upstream.txt; cat Upstream.txt | sed 's/.../& /2g;s/./& /g' > Upstream.tab
	grep -v '>' Downstream.fa > Downstream.txt; cat Downstream.txt | sed 's/.../& /2g;s/./& /g' > Downstream.tab

	for i in {1..100}; do
		#Location of output files
		UpstreamLists=$sample.Column.$i.Upstream.$reference.$subset.txt
		DownstreamLists=$sample.Column.$i.Downstream.$reference.$subset.txt
		
		#Save lists of dNMPs at each upstream/downstream position
		awk -v field=$i '{ print $field }' Upstream.tab > $UpstreamLists
		awk -v field=$i '{ print $field }' Downstream.tab > $DownstreamLists
	done
	
#############################################################################################################################
	#STEP 5: Calculate frequencies of dNMPs +/- 100 base pairs from rNMPs

	#Calculate frequencies at each position
	for direction in "Upstream" "Downstream"; do
		
		for file in `ls -v ./$sample.Column.*.$direction.$reference.$subset.txt`; do
		
		#Calculate count of each dNMP
		A_FlankCount=$(grep -o 'A' $file | wc -l); C_FlankCount=$(grep -o 'C' $file | wc -l)
		G_FlankCount=$(grep -o 'G' $file | wc -l); T_FlankCount=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNMPs
		FlankCount=$(($A_FlankCount+$C_FlankCount+$G_FlankCount+$T_FlankCount))

		#Calculate normalized frequencies of dNMPs
		A_FlankFreq=$(echo "scale=12; ($A_FlankCount/$FlankCount)/$A_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		C_FlankFreq=$(echo "scale=12; ($C_FlankCount/$FlankCount)/$C_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		G_FlankFreq=$(echo "scale=12; ($G_FlankCount/$FlankCount)/$G_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		T_FlankFreq=$(echo "scale=12; ($T_FlankCount/$FlankCount)/$T_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNMPs frequencies to TXT file
		echo $A_FlankFreq >> A_FlankFreq.$direction.txt; echo $C_FlankFreq >> C_FlankFreq.$direction.txt
		echo $G_FlankFreq >> G_FlankFreq.$direction.txt; echo $T_FlankFreq >> T_FlankFreq.$direction.txt
		
		if [ $i == "Downstream" ]; then
			DownFreq=$(paste A_FlankFreq.Downstream.txt C_FlankFreq.Downstream.txt G_FlankFreq.Downstream.txt T_FlankFreq.Downstream.txt)
		elif [ $i == "Upstream" ]; then
			UpFreq=$(paste A_FlankFreq.Upstream.txt C_FlankFreq.Upstream.txt G_FlankFreq.Upstream.txt T_FlankFreq.Upstream.txt | tac -)
		fi
				
		done
	done
	
#############################################################################################################################
	#STEP 6: Create and save dataset file containing nucleotide frequencies

	#Add nucleotide to header line
	echo -e "\tA\tC\tG\tU/T" > $dataset
	
	#Add nucleotide positions and frequencies in correct order
	Freqs=$(cat <(echo "$UpFreq") <(echo "$RiboFreq") <(echo "$DownFreq"))
	paste <(echo "$(seq -100 1 100)") <(cat <(echo "$Freqs")) >> $dataset

	#Remove temp files
	rm -f ./*Freq.txt ./*Upstream.* ./*Downstream.* ./RiboBases.txt ./temp.fa*

	#Let the user know the analysis is complete
	echo "Calculation of nucleotide frequencies for $sample ($subset) is complete"
	
	done
done
