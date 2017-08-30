#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates rNMP frequencies and flanking dNMP frequencies (+/- 100 bp)

#Usage statement
function usage () {
	echo "Usage: Frequencies.sh [options]
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

for sample in ${sample[@]}; do
	for subset in "mito" "nucleus"; do
	
#############################################################################################################################
	#Input files
	reads=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-ReadInformation.$subset.txt
	BED=$directory/Ribose-Map/References/$reference.bed; FASTA=$directory/Ribose-Map/References/$reference.fa
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$sample-Coordinates.bed
	
	#Output directory
	output=$directory/Ribose-Map/Results/$reference/$sample/Frequencies
	
	#Create directory
	mkdir -p $output
	
	if [ -s $coordinates ]; then

	#Remove old file
	rm -f $output/$sample-*.txt
	
#############################################################################################################################
	#STEP 1: Calculate frequencies of reference genome
	
	#Subset FASTA file based on region
	if [ $subset == "mito" ]; then
		chr=$(awk '{print $1}' $BED | grep -E '(chrM|MT)')
		samtools faidx $FASTA $chr > $output/temp.fa; samtools faidx $output/temp.fa
	elif [ $subset == "nucleus" ]; then
		chr=$(awk '{print $1}' $BED | grep -vE '(chrM|MT)')
		samtools faidx $FASTA $chr > $output/temp.fa; samtools faidx $output/temp.fa
	fi

	#Calculate counts of each nucleotide
	A_BkgCount=$(grep -v '>' $output/temp.fa | grep -o 'A' - | wc -l)
	C_BkgCount=$(grep -v '>' $output/temp.fa | grep -o 'C' - | wc -l)
	G_BkgCount=$(grep -v '>' $output/temp.fa | grep -o 'G' - | wc -l)
	T_BkgCount=$(grep -v '>' $output/temp.fa | grep -o 'T' - | wc -l)
	
	#Calculate total number of nucleotides
	total_Bkg=$(($A_BkgCount+$C_BkgCount+$G_BkgCount+$T_BkgCount))
	
	#Calculate frequency of each nucleotide
	A_BkgFreq=$(echo "scale=12; ($A_BkgCount+$T_BkgCount)/($total_Bkg*2)" | bc | awk '{printf "%.12f\n", $0}')
	C_BkgFreq=$(echo "scale=12; ($C_BkgCount+$G_BkgCount)/($total_Bkg*2)" | bc | awk '{printf "%.12f\n", $0}')
	G_BkgFreq=$(echo "scale=12; ($G_BkgCount+$C_BkgCount)/($total_Bkg*2)" | bc | awk '{printf "%.12f\n", $0}')
	T_BkgFreq=$(echo "scale=12; ($T_BkgCount+$A_BkgCount)/($total_Bkg*2)" | bc | awk '{printf "%.12f\n", $0}')
	
#############################################################################################################################
	#STEP 2: Calculate frequencies of rNMPs in libraries
	
	#Subset and sort coordinates based on genomic region
	if [ $subset == "mito" ]; then
		grep -E '(chrM|MT)' $coordinates > $output/$sample-Coordinates.$subset.bed
	elif [ $subset == "nucleus" ]; then
		grep -vE '(chrM|MT)' $coordinates > $output/$sample-Coordinates.$subset.bed
	fi
	
	#Extract rNMP bases
	bedtools getfasta -s -fi $output/temp.fa -bed $output/$sample-Coordinates.$subset.bed \
	| grep -v '>' - > $output/RiboBases.txt
	
	#Calculate counts of rNMPs
	A_RiboCount=$(awk '$1 == "A"' $output/RiboBases.txt | wc -l)
	C_RiboCount=$(awk '$1 == "C"' $output/RiboBases.txt | wc -l)
	G_RiboCount=$(awk '$1 == "G"' $output/RiboBases.txt | wc -l)
	U_RiboCount=$(awk '$1 == "T"' $output/RiboBases.txt | wc -l)
	
	#Calculate total number of rNMPs
	RiboCount=$(($A_RiboCount+$C_RiboCount+$G_RiboCount+$U_RiboCount))

	#Calculate normalized frequency of each rNMP
	A_RiboFreq=$(echo "scale=12; ($A_RiboCount/$RiboCount)/$A_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	C_RiboFreq=$(echo "scale=12; ($C_RiboCount/$RiboCount)/$C_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	G_RiboFreq=$(echo "scale=12; ($G_RiboCount/$RiboCount)/$G_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	U_RiboFreq=$(echo "scale=12; ($U_RiboCount/$RiboCount)/$T_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')

	#Save normalized frequencies of rNMPs together
	Ribo=$(echo -e "$A_RiboFreq\t$C_RiboFreq\t$G_RiboFreq\t$U_RiboFreq")
	
#############################################################################################################################
	#STEP 3: Obtain coordinates/sequences of dNMPs +/- 100 bp from rNMPs

	#Obtain coordinates of flanking sequences and remove coordinates where start = end
	bedtools flank -i $coordinates -s -g $BED -l 100 -r 0 | awk '$2 != $3' - > $output/Up.bed
	bedtools flank -i $coordinates -s -g $BED -l 0 -r 100 | awk '$2 != $3' - > $output/Down.bed
	
	#Obtain nucleotide sequences flanking rNMPs using coordinates from above
	bedtools getfasta -s -fi temp.fa -bed $output/Upstream.bed -fo $output/Up.fa
	bedtools getfasta -s -fi temp.fa -bed $output/Downstream.bed -fo $output/Down.fa
	
#############################################################################################################################
	#STEP 4: Insert tabs between sequences of dNMPs +/- 100 bp from rNMPs

	#Extract sequences (reverse order of upstream)
	grep -v '>' $output/Upstream.fa | rev > $output/Up.txt
	grep -v '>' $output/Downstream.fa > $output/Down.txt
	
	#Insert tabs between each base for easier parsing
	cat $output/Upstream.txt | sed 's/.../& /2g;s/./& /g' > $output/Up.tab
	cat $output/Downstream.txt | sed 's/.../& /2g;s/./& /g' > $output/Down.tab

	#Save lists of dNMPs at each of the +/-100 positions in separate files
	for i in {1..100}; do
		awk -v field=$i '{ print $field }' $output/Up.tab > $output/$sample.Up.$i.txt
		awk -v field=$i '{ print $field }' $output/Down.tab > $output/$sample.Down.$i.txt
	done
	
#############################################################################################################################
	#STEP 5: Calculate frequencies of dNMPs +/- 100 base pairs from rNMPs

	for direction in "Up" "Down"; do
		
		#'-v' = natural sort of #'s
		for file in `ls -v $output/$sample.$direction.$i.txt`; do
		
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
		echo $A_FlankFreq >> $output/A_$direction.txt; echo $C_FlankFreq >> $output/C_$direction.txt
		echo $G_FlankFreq >> $output/G_$direction.txt; echo $T_FlankFreq >> $output/T_$direction.txt
		
		if [ $direction == "Up" ]; then
			Up=$(paste $output/A_Up.txt $output/C_Up.txt $output/G_Up.txt $output/T_Up.txt | tac -)
		elif [ $direction == "Down" ]; then
			Down=$(paste $output/A_Down.txt $output/C_Down.txt $output/G_Down.txt $output/T_Down.txt)
		fi
				
		done
	done
	
#############################################################################################################################
	#STEP 6: Create and save dataset file containing nucleotide frequencies

	#Add nucleotides to header line
	echo -e "\tA\tC\tG\tU/T" > $output/$sample-Frequencies.$subset.txt

	#Add positions and frequencies of nucleotides in correct order to create dataset
	paste <(echo "$(seq -100 1 100)") <(cat <(echo "$Up") <(echo "$Ribo") <(echo "$Down")) \
	>> $output/$sample-Frequencies.$subset.txt

#############################################################################################################################
	#STEP 7: Create and save file containing background nucleotide frequencies
		
	#Add nucleotides to header line
	echo -e "A\tC\tG\tT" > $output/$sample-BackgroundFrequencies.txt
	
	#Add frequencies of nucleotides in reference genome
	paste <(echo -e "$A_BkgFreq\t$C_BkgFreq\t$G_BkgFreq\t$T_BkgFreq") >> $output/$sample-BackgroundFrequencies.txt
	
	#Add total number of nucleotides in reference genome
	echo -e "Total number of bases in reference genome: $((total_Bkg*2))" >> $output/$sample-BackgroundFrequencies.txt

#############################################################################################################################
	#Print completion status
	echo "Calculation of frequencies for $sample ($subset) is complete"
	
	#Remove temp files
	rm -f $output/*Up.* $output/*Down.* $output/RiboBases.txt $output/temp.fa*
	
	fi
	
	done
done
