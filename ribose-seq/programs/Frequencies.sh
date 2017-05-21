#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates rNMP frequencies and flanking dNMP frequencies (+/- 100 bp)

#Usage statement
function usage () {
	echo "Usage: Frequencies.sh [-i] 'Sample(s)' [-r] 'Reference' [-d] 'Directory' [-h]
	-i Input sample(s) (e.g., FS1, FS2, FS3 etc.)
	-r Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)
	-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "i:r:d:h" opt; do
    case $opt in
        #Allow multiple input arguments
        i ) sample=($OPTARG) ;;
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
	for subset in "all" "mito" "nucleus"; do

#############################################################################################################################
	#Input files
	reads=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$subset/$sample-ReadInformation.$subset.txt
	BED=$directory/Ribose-Map/Reference/$reference.bed; FASTA=$directory/Ribose-Map/Reference/$reference.$subset.fa
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$subset/$sample-Coordinates.$subset.bed

	#Output directory and file
	output=$directory/Ribose-Map/Results/$reference/$sample/Frequencies/$subset
	dataset=$output/$sample-Frequencies.$reference.$subset.txt
	
	#Create directory and remove old file
	mkdir -p $output; rm -f $output/{*.txt}
		
#############################################################################################################################
	#STEP 1: Calculate frequencies of reference genome
	
	#Calculate counts of each nucleotide
	A_BkgCount=$(grep -v '>' $FASTA | grep -o 'A' - | wc -l); C_BkgCount=$(grep -v '>' $FASTA | grep -o 'C' - | wc -l)
	G_BkgCount=$(grep -v '>' $FASTA | grep -o 'G' - | wc -l); T_BkgCount=$(grep -v '>' $FASTA | grep -o 'T' - | wc -l)
	
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
		cat $reads | awk '{print substr($0,length($0))}' - > riboSequences.txt
	elif [ $subset == "mito" ]; then
		#Select only reads located in mito DNA and extract rNMP bases
		grep -E '(chrM|MT)' $reads | awk '{print substr($0,length($0))}' - > riboSequences.txt
	elif [ $subset == "nucleus" ]; then
		#Select only reads located in nuclear DNA and extract rNMP bases
		grep -v -E '(chrM|MT)' $reads | awk '{print substr($0,length($0))}' - > riboSequences.txt
	fi
	
	#Calculate counts of rNMPs
	A_RiboCount=$(awk '$1 == "A"' riboSequences.txt | wc -l); C_RiboCount=$(awk '$1 == "C"' riboSequences.txt | wc -l)
	G_RiboCount=$(awk '$1 == "G"' riboSequences.txt | wc -l); U_RiboCount=$(awk '$1 == "T"' riboSequences.txt | wc -l)
	
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
	bedtools getfasta -s -fi $FASTA -bed Upstream.bed -fo Upstream.fasta
	bedtools getfasta -s -fi $FASTA -bed Downstream.bed -fo Downstream.fasta

#############################################################################################################################
	#STEP 4: Insert tabs between sequences of dNMPs +/- 100 bp from rNMPs

	#Extract sequences from FASTA files (reverse order of upstream) and insert tabs between each nucleotide
	grep -v '>' Upstream.fasta | rev > Upstream.txt; cat Upstream.txt | sed 's/.../& /2g;s/./& /g' > Upstream.tab
	grep -v '>' Downstream.fasta > Downstream.txt; cat Downstream.txt | sed 's/.../& /2g;s/./& /g' > Downstream.tab

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
	for file in `ls -v ./$sample.Column.*.Upstream.$reference.$subset.txt`; do

		#Calculate count of each dNMP
		A_UpCount=$(grep -o 'A' $file | wc -l); C_UpCount=$(grep -o 'C' $file | wc -l)
		G_UpCount=$(grep -o 'G' $file | wc -l); T_UpCount=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNMPs
		UpCount=$(($A_UpCount+$C_UpCount+$G_UpCount+$T_UpCount))

		#Calculate normalized frequencies of dNMPs
		A_UpFreq=$(echo "scale=12; ($A_UpCount/$UpCount)/$A_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		C_UpFreq=$(echo "scale=12; ($C_UpCount/$UpCount)/$C_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		G_UpFreq=$(echo "scale=12; ($G_UpCount/$UpCount)/$G_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		T_UpFreq=$(echo "scale=12; ($T_UpCount/$UpCount)/$T_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNMPs frequencies to TXT file
		echo $A_UpFreq >> A_UpFreq.txt; echo $C_UpFreq >> C_UpFreq.txt
		echo $G_UpFreq >> G_UpFreq.txt; echo $T_UpFreq >> T_UpFreq.txt
			
		#Combine upstream dNMP frequencies and reverse (order = -100 --> -1)
		UpFreq=$(paste A_UpFreq.txt C_UpFreq.txt G_UpFreq.txt T_UpFreq.txt | tac -)
		
	done

	#Calculate frequencies at each position
	for file in `ls -v ./$sample.Column.*.Downstream.$reference.$subset.txt`; do

		#Calculate count of each dNMP
		A_DownCount=$(grep -o 'A' $file | wc -l); C_DownCount=$(grep -o 'C' $file | wc -l)
		G_DownCount=$(grep -o 'G' $file | wc -l); T_DownCount=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNMPs
		DownCount=$(($A_DownCount+$C_DownCount+$G_DownCount+$T_DownCount))
	
		#Calculate normalized frequencies of dNMPs
		A_DownFreq=$(echo "scale=12; ($A_DownCount/$DownCount)/$A_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		C_DownFreq=$(echo "scale=12; ($C_DownCount/$DownCount)/$C_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		G_DownFreq=$(echo "scale=12; ($G_DownCount/$DownCount)/$G_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		T_DownFreq=$(echo "scale=12; ($T_DownCount/$DownCount)/$T_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNMPs frequencies to TXT file
		echo $A_DownFreq >> A_DownFreq.txt; echo $C_DownFreq >> C_DownFreq.txt
		echo $G_DownFreq >> G_DownFreq.txt; echo $T_DownFreq >> T_DownFreq.txt
		
		#Combine downstream dNMP frequencies (order = +1 --> +100)
		DownFreq=$(paste A_DownFreq.txt C_DownFreq.txt G_DownFreq.txt T_DownFreq.txt)
	
	done
	
#############################################################################################################################
	#STEP 6: Create and save dataset file containing nucleotide frequencies

	#Add nucleotide to header line
	echo -e "\tA\tC\tG\tU/T" > $dataset
	
	#Add nucleotide positions and frequencies in correct order
	Freqs=$(cat <(echo "$UpFreq") <(echo "$RiboFreq") <(echo "$DownFreq"))
	paste <(echo "$(seq -100 1 100)") <(cat <(echo "$Freqs")) >> $dataset

	#Let the user know the analysis is complete
	echo "Calculation of nucleotide frequencies for $sample ($subset) is complete"

done
done

#Remove temp files
rm -f ./*Freq.txt ./*Upstream.* ./*Downstream.* ./riboSequences.txt
