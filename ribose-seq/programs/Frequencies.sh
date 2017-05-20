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

subset=("all" "mito" "nucleus")

#Calculate frequencies
for sample in ${sample[@]}; do
	for subset in ${subset[@]}; do

#############################################################################################################################
	#Input files
	reads=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$subset/$sample-ReadInformation.$subset.txt
	BED=$directory/Ribose-Map/Reference/$reference.bed; FASTA=$directory/Ribose-Map/Reference/$reference.$subset.fa
	coordinates=$directory/Ribose-Map/Results/$reference/$sample/Coordinates/$subset/$sample-Coordinates.$subset.bed

	#Output directories
	output1=$directory/Ribose-Map/Results/$reference/$sample/Frequencies/Datasets/$subset
	output2=$directory/Ribose-Map/Results/$reference/$sample/Frequencies/dNMPs/$subset/Columns/upstream
	output3=$directory/Ribose-Map/Results/$reference/$sample/Frequencies/dNMPs/$subset/Columns/downstream

	#Create directories and remove older versions files
	mkdir -p $output{1..3}; rm -f $output{1..3}/*.txt $output{1..3}/*.tab
	
	#Output files
	dataset1=$output1/$sample-NucleotideFrequenciesv1.$reference.$subset.txt
	dataset2=$output1/$sample-NucleotideFrequenciesv2.$reference.$subset.txt
		
#############################################################################################################################
	#STEP 1: Calculate background nucleotide frequencies of reference genome
	
	#Calculate counts of each nucleotide
	A_BkgCount=$(grep -v '>' $FASTA | grep -o 'A' - | wc -l)
	C_BkgCount=$(grep -v '>' $FASTA | grep -o 'C' - | wc -l)
	G_BkgCount=$(grep -v '>' $FASTA | grep -o 'G' - | wc -l)
	T_BkgCount=$(grep -v '>' $FASTA | grep -o 'T' - | wc -l)
	
	#Calculate total number of nucleotides
	total_Bkg=$(($A_BkgCount+$C_BkgCount+$G_BkgCount+$T_BkgCount))

	#Calculate frequency of each nucleotide
	A_BkgFreq=$(echo "scale = 12; $A_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')
	C_BkgFreq=$(echo "scale = 12; $C_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')
	G_BkgFreq=$(echo "scale = 12; $G_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')
	T_BkgFreq=$(echo "scale = 12; $T_BkgCount/$total_Bkg" | bc | awk '{printf "%.12f\n", $0}')

#############################################################################################################################
	#STEP 2: Calculate rNMP Frequencies

	if [ $subset == "all" ]; then
		#Select all reads located in genomic DNA and extract rNMP sequences
		cat $reads | awk '{print substr($0,length($0))}' - > riboSequences.txt
	elif [ $subset == "mito" ]; then
		#Select only reads located in mitochondrial DNA and extract rNMP sequences
		grep -E '(chrM|MT)' $reads | awk '{print substr($0,length($0))}' - > riboSequences.txt
	elif [ $subset == "nucleus" ]; then
		#Select only reads located in nuclear DNA and extract rNMP sequences
		grep -v -E '(chrM|MT)' $reads | awk '{print substr($0,length($0))}' - > riboSequences.txt
	fi
	
	#Calculate counts of rNMPs
	A_Count1=$(awk '$1 == "A"' riboSequences.txt | wc -l)
	C_Count1=$(awk '$1 == "C"' riboSequences.txt | wc -l)
	G_Count1=$(awk '$1 == "G"' riboSequences.txt | wc -l)
	U_Count1=$(awk '$1 == "T"' riboSequences.txt | wc -l)
	
	#Calculate total number of rNMPs
	total1=$(($A_Count1+$C_Count1+$G_Count1+$U_Count1))

	#Calculate normalized frequency of each rNMP
	A_Frequency1=$(echo "scale = 12; ($A_Count1/$total1)/$A_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	C_Frequency1=$(echo "scale = 12; ($C_Count1/$total1)/$C_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	G_Frequency1=$(echo "scale = 12; ($G_Count1/$total1)/$G_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')
	U_Frequency1=$(echo "scale = 12; ($U_Count1/$total1)/$T_BkgFreq" | bc | awk '{printf "%.12f\n", $0}')

	#Save normalized frequencies of rNMPs together
	riboFrequencies=$(echo -e "$A_Frequency1\t$C_Frequency1\t$G_Frequency1\t$U_Frequency1")
	
#############################################################################################################################
	#STEP 3: Obtain coordinates and sequences of +/- 100 downstream/upstream dNMPs from rNMPs

	#Obtain coordinates of sequences upstream/downstream from rNMPs
	bedtools flank -i $coordinates -s -g $BED -l 100 -r 0 | awk '$2 != "0" && $3 != "0"' - > upstreamIntervals.txt
	bedtools flank -i $coordinates -s -g $BED -l 0 -r 100 | awk '$2 != "0" && $3 != "0"' - > downstreamIntervals.txt
	
	#Obtain sequences of the sequences upstream/downstream from rNMPs
	bedtools getfasta -s -fi $FASTA -bed upstreamIntervals.txt -fo upstreamSequences.txt
	bedtools getfasta -s -fi $FASTA -bed downstreamIntervals.txt -fo downstreamSequences.txt

#############################################################################################################################
	#STEP 4: Tabulate sequences of dNMPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Extract sequences from FASTA files
	#Reverse order of upstream nucleotides
	grep -v '>' upstreamSequences.txt | rev > sequences1.txt
	grep -v '>' downstreamSequences.txt > sequences2.txt
				
	#Insert tabs between each nucleotide
	cat sequences1.txt | sed 's/.../& /2g;s/./& /g' > columns1.tab
	cat sequences2.txt | sed 's/.../& /2g;s/./& /g' > columns2.tab

	for i in {1..100}; do
		#Location of output files
		lists1=$output3/$sample.column.$i.upstream.$reference.$subset.txt
		lists2=$output4/$sample.column.$i.downstream.$reference.$subset.txt
		#Save lists of dNMPs at each +/- 100 bp downstream/upstream position
		awk -v field=$i '{ print $field }' columns1.tab > $lists1
		awk -v field=$i '{ print $field }' columns2.tab > $lists2
	done

#############################################################################################################################
	#STEP 5: Calculate frequencies of dNMPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Calculate frequencies at each position
	for file in `ls -v $output3/$sample*.txt`; do

		#Calculate count of each dNMP
		A_Count2=$(grep -o 'A' $file | wc -l)
		C_Count2=$(grep -o 'C' $file | wc -l)
		G_Count2=$(grep -o 'G' $file | wc -l)
		T_Count2=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNMPs
		total2=$(($A_Count2+$C_Count2+$G_Count2+$T_Count2))

		#Calculate normalized frequencies of dNMPs
		A_Frequency2=$(echo "scale = 12; ($A_Count2/$total2)/$A_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		C_Frequency2=$(echo "scale = 12; ($C_Count2/$total2)/$C_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		G_Frequency2=$(echo "scale = 12; ($G_Count2/$total2)/$G_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		T_Frequency2=$(echo "scale = 12; ($T_Count2/$total2)/$T_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNMPs frequencies to TXT file
		echo $A_Frequency2 >> A_frequency2.txt; echo $C_Frequency2 >> C_frequency2.txt
		echo $G_Frequency2 >> G_frequency2.txt; echo $T_Frequency2 >> T_frequency2.txt
			
		#Combine upstream dNMP frequencies together and reverse (frequencies ordered from -100 --> -1)
		upstreamFrequencies=$(paste A_frequency2.txt C_frequency2.txt G_frequency2.txt T_frequency2.txt | tac -)
		
	done

	#Calculate frequencies at each position
	for file in `ls -v $output4/$sample*.txt`; do

		#Calculate count of each dNMP
		A_Count3=$(grep -v '>' $file | grep -o 'A' - | wc -l)
		C_Count3=$(grep -v '>' $file | grep -o 'C' - | wc -l)
		G_Count3=$(grep -v '>' $file | grep -o 'G' - | wc -l)
		T_Count3=$(grep -v '>' $file | grep -o 'T' - | wc -l)

		#Calculate total number of dNMPs
		total3=$(($A_Count3+$C_Count3+$G_Count3+$T_Count3))
	
		#Calculate normalized frequencies of dNMPs
		A_Frequency3=$(echo "scale = 12; ($A_Count3/$total3)/$A_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		C_Frequency3=$(echo "scale = 12; ($C_Count3/$total3)/$C_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		G_Frequency3=$(echo "scale = 12; ($G_Count3/$total3)/$G_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		T_Frequency3=$(echo "scale = 12; ($T_Count3/$total3)/$T_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNMPs frequencies to TXT file
		echo $A_Frequency3 >> A_frequency3.txt; echo $C_Frequency3 >> C_frequency3.txt
		echo $G_Frequency3 >> G_frequency3.txt; echo $T_Frequency3 >> T_frequency3.txt
		
		#Combine downstream dNMP frequencies together as is (frequencies ordered from +1 --> +100)
		downstreamFrequencies=$(paste A_frequency3.txt C_frequency3.txt G_frequency3.txt T_frequency3.txt)
	
	done

	#Remove intermediate files
	rm A_frequency{2..3}.txt C_frequency{2..3}.txt G_frequency{2..3}.txt T_frequency{2..3}.txt
	
#############################################################################################################################
	#STEP 6: Create dataset file containing nucleotide frequencies for plotting

	#Combine rNMP frequencies and upstream and downstream dNMP frequencies in appropriate order
	data1=$(cat <(echo "$upstreamFrequencies") <(echo "$riboFrequencies") <(echo "$downstreamFrequencies"))
	
	#Add nucleotide positions (-100 --> +100) and nucleotide symbols to header line (A, C, G, and U/T)
	echo -e "\tA\tC\tG\tU/T" > $dataset1 && paste <(echo "$(seq -100 1 100)") <(cat <(echo "$data1")) >> $dataset1

	#Smaller dataset (-15 nt to +15 nt)
	data2=$(head -117 $dataset1 | tail -31)
	echo -e "\tA\tC\tG\tU/T" > $dataset2 && cat <(echo "$data2") >> $dataset2

	#Let the user know the program is has finished running
	echo "Calculation of nucleotide frequencies for $sample ($subset) is complete"

done
done

#Remove temp files
rm -f ./*.txt ./*.tab
