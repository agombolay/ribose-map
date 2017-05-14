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
	-d Directory (e.g., /projects/home/agombolay3/data/repository/Ribose-seq-Project)"
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

#subset=("genome" "nucleus" "mitochondria")
subset=("mitochondria")

#Calculate frequencies
for sample in ${sample[@]}; do
	for subset in ${subset[@]}; do

#############################################################################################################################
	#Input files
	referenceBED=$directory/ribose-seq/reference/$reference.bed
	referenceFasta1=$directory/ribose-seq/reference/$reference.fa
	referenceFasta2=$directory/ribose-seq/reference/$reference.$subset.fa
	
	reads=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset/$sample.read-information.$subset.txt
	coordinates=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset/$sample.rNMP-coordinates.$subset.bed

	#Output directories
	output1=$directory/ribose-seq/results/Background-Frequencies
	output2=$directory/ribose-seq/results/$reference/$sample/Frequencies/Datasets/$subset
	output3=$directory/ribose-seq/results/$reference/$sample/Frequencies/dNMPs/$subset/Columns/upstream
	output4=$directory/ribose-seq/results/$reference/$sample/Frequencies/dNMPs/$subset/Columns/downstream

	#Create directories
	mkdir -p $output{1..4}

	#Remove older versions files
	rm -f $output{1..4}/*.txt $output{1..4}/*.tab
	
	#Output files
	background=$output1/Background-Frequencies.$reference.$subset.txt
	zoomed=$output2/$sample.nucleotide-frequencies-zoomed.$reference.$subset.txt
	dataset=$output2/$sample.nucleotide-frequencies-dataset.$reference.$subset.txt
		
#############################################################################################################################
	#STEP 1: Calculate background dNMP frequencies of reference genome

	#Index reference FASTA file
	samtools faidx $referenceFasta1
	
	if [ ! -f $referenceFasta2 ]; then

		#Specify all genomic DNA
		if [ $subset == "genome" ]; then
			cp $referenceFasta1 $referenceFasta2
		
		#Subset mitochondrial DNA
		elif [ $reference == "pombe" ] && [ $subset == "mitochondria" ]; then
			samtools faidx $referenceFasta1 MT > $referenceFasta2
		
		elif [ $reference != "pombe" ] && [ $subset == "mitochondria" ]; then
			samtools faidx $referenceFasta1 chrM > $referenceFasta2
			
		#Subset nuclear DNA
		elif [ $reference == "pombe" ] && [ $subset == "nucleus" ]; then
			samtools faidx $referenceFasta1 I II III > $referenceFasta2
		
		elif [ $reference == "sacCer2" ] && [ $subset == "nucleus" ]; then
			samtools faidx $referenceFasta1 2micron chrI chrII chrIII chrIV chrV chrVI chrVII \
			chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI > $referenceFasta2
		
		elif [ $reference == "mm9" ] && [ $subset == "nucleus" ]; then chr $(seq 1 1 19)" X Y";
			for i in $chr; do samtools faidx $referenceFasta1 chr$i > $referenceFasta2; done
		
		elif [ $reference == "hg38" ] && [ $subset == "nucleus" ]; then chr $(seq 1 1 22)" X Y";
			for i in $chr; do samtools faidx $referenceFasta1 chr$i > $referenceFasta2; done
		fi
	fi
	
	#Index reference FASTA file
	samtools faidx $referenceFasta2
	
	#Calculate counts of each dNMP
	A_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'A' - | wc -l)
	C_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'C' - | wc -l)
	G_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'G' - | wc -l)
	T_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'T' - | wc -l)
	
	#Calculate total number of dNMPs
	total0=$(($A_Count0+$C_Count0+$G_Count0+$T_Count0))

	#Calculate frequency of each dNMP
	A_Frequency0=$(echo "scale = 12; $A_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')
	C_Frequency0=$(echo "scale = 12; $C_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')
	G_Frequency0=$(echo "scale = 12; $G_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')
	T_Frequency0=$(echo "scale = 12; $T_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')

	#Save frequencies of dNMPs in the reference FASTA file (background frequencies) to TXT file
	echo -e "A\tC\tG\tU/T\n$A_Frequency0\t$C_Frequency0\t$G_Frequency0\t$T_Frequency0" > $background

#############################################################################################################################
	#STEP 2: Calculate rNMP Frequencies

	if [ $subset == "genome" ]; then
		#Select all reads located in genomic DNA and extract rNMP sequences
		cat $reads | awk '{print substr($0,length($0))}' - > riboSequences.txt
	elif [ $subset == "mitochondria" ]; then
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
	A_Frequency1=$(echo "scale = 12; ($A_Count1/$total1)/$A_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	C_Frequency1=$(echo "scale = 12; ($C_Count1/$total1)/$C_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	G_Frequency1=$(echo "scale = 12; ($G_Count1/$total1)/$G_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	U_Frequency1=$(echo "scale = 12; ($U_Count1/$total1)/$T_Frequency0" | bc | awk '{printf "%.12f\n", $0}')

	#Save normalized frequencies of rNMPs together
	riboFrequencies=$(echo -e "$A_Frequency1\t$C_Frequency1\t$G_Frequency1\t$U_Frequency1")
	
#############################################################################################################################
	#STEP 3: Obtain coordinates and sequences of +/- 100 downstream/upstream dNMPs from rNMPs

	#Obtain coordinates of upstream/downstream sequences based on rNMP coordinates
	#bedtools flank -i $coordinates -s -g $referenceBED -l 100 -r 0 > upstreamIntervals.txt
	#bedtools flank -i $coordinates -s -g $referenceBED -l 0 -r 100 > downstreamIntervals.txt
	
	#Obtain sequences of upstream/downstream coordinates
	#bedtools getfasta -s -fi $referenceFasta1 -bed upstreamIntervals.txt -fo upstreamSequences.txt
	#bedtools getfasta -s -fi $referenceFasta1 -bed downstreamIntervals.txt -fo downstreamSequences.txt

#############################################################################################################################
	#STEP 4: Tabulate sequences of dNMPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Extract sequences from FASTA files
	#Reverse order of upstream nucleotides
	#grep -v '>' upstreamSequences.txt | rev > sequences1.txt
	#grep -v '>' downstreamSequences.txt > sequences2.txt
				
	#Insert tabs between each nucleotide
	#cat sequences1.txt | sed 's/.../& /2g;s/./& /g' > columns1.tab
	#cat sequences2.txt | sed 's/.../& /2g;s/./& /g' > columns2.tab

	#for i in {1..100}; do
		#Location of output files
	#	lists1=$output3/$sample.column.$i.upstream.$reference.$subset.txt
	#	lists2=$output4/$sample.column.$i.downstream.$reference.$subset.txt
		#Save lists of dNMPs at each +/- 100 bp downstream/upstream position
	#	awk -v field=$i '{ print $field }' columns1.tab > $lists1
	#	awk -v field=$i '{ print $field }' columns2.tab > $lists2
	#done

#############################################################################################################################
	#STEP 5: Calculate frequencies of dNMPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Calculate frequencies at each position
	#for file in `ls -v $output3/$sample*.txt`; do

		#Calculate count of each dNMP
	#	A_Count2=$(grep -o 'A' $file | wc -l)
	#	C_Count2=$(grep -o 'C' $file | wc -l)
	#	G_Count2=$(grep -o 'G' $file | wc -l)
	#	T_Count2=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNMPs
	#	total2=$(($A_Count2+$C_Count2+$G_Count2+$T_Count2))

		#Calculate normalized frequencies of dNMPs
	#	A_Frequency2=$(echo "scale = 12; ($A_Count2/$total2)/$A_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	#	C_Frequency2=$(echo "scale = 12; ($C_Count2/$total2)/$C_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	#	G_Frequency2=$(echo "scale = 12; ($G_Count2/$total2)/$G_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	#	T_Frequency2=$(echo "scale = 12; ($T_Count2/$total2)/$T_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNMPs frequencies to TXT file
	#	echo $A_Frequency2 >> A_frequency2.txt; echo $C_Frequency2 >> C_frequency2.txt
	#	echo $G_Frequency2 >> G_frequency2.txt; echo $T_Frequency2 >> T_frequency2.txt
			
		#Combine upstream dNMP frequencies together and reverse (frequencies ordered from -100 --> -1)
	#	upstreamFrequencies=$(paste A_frequency2.txt C_frequency2.txt G_frequency2.txt T_frequency2.txt | tac -)
		
	#done

	#Calculate frequencies at each position
	#for file in `ls -v $output4/$sample*.txt`; do

		#Calculate count of each dNMP
	#	A_Count3=$(grep -v '>' $file | grep -o 'A' - | wc -l)
	#	C_Count3=$(grep -v '>' $file | grep -o 'C' - | wc -l)
	#	G_Count3=$(grep -v '>' $file | grep -o 'G' - | wc -l)
	#	T_Count3=$(grep -v '>' $file | grep -o 'T' - | wc -l)

		#Calculate total number of dNMPs
	#	total3=$(($A_Count3+$C_Count3+$G_Count3+$T_Count3))
	
		#Calculate normalized frequencies of dNMPs
	#	A_Frequency3=$(echo "scale = 12; ($A_Count3/$total3)/$A_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	#	C_Frequency3=$(echo "scale = 12; ($C_Count3/$total3)/$C_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	#	G_Frequency3=$(echo "scale = 12; ($G_Count3/$total3)/$G_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
	#	T_Frequency3=$(echo "scale = 12; ($T_Count3/$total3)/$T_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNMPs frequencies to TXT file
	#	echo $A_Frequency3 >> A_frequency3.txt; echo $C_Frequency3 >> C_frequency3.txt
	#	echo $G_Frequency3 >> G_frequency3.txt; echo $T_Frequency3 >> T_frequency3.txt
		
		#Combine downstream dNMP frequencies together as is (frequencies ordered from +1 --> +100)
	#	downstreamFrequencies=$(paste A_frequency3.txt C_frequency3.txt G_frequency3.txt T_frequency3.txt)
	
	#done

	#Remove intermediate files
	#rm A_frequency{2..3}.txt C_frequency{2..3}.txt G_frequency{2..3}.txt T_frequency{2..3}.txt
	
#############################################################################################################################
	#STEP 6: Create dataset file containing nucleotide frequencies for plotting

	#Combine rNMP frequencies and upstream and downstream dNMP frequencies in appropriate order
	#data1=$(cat <(echo "$upstreamFrequencies") <(echo "$riboFrequencies") <(echo "$downstreamFrequencies"))
	
	#Add nucleotide positions (-100 --> +100) and nucleotide symbols to header line (A, C, G, and U/T)
	#echo -e "\tA\tC\tG\tU/T" > $dataset && paste <(echo "$(seq -100 1 100)") <(cat <(echo "$data1")) >> $dataset

	#Smaller dataset (-15 nt to +15 nt)
	#data2=$(head -117 $dataset | tail -31)
	#echo -e "\tA\tC\tG\tU/T" > $zoomed && cat <(echo "$data2") >> $zoomed

	#Let the user know the program is has finished running
	#echo "Calculation of nucleotide frequencies for $sample ($subset) is complete"

done
done

#Remove temp files
#rm -f ./*.txt ./*.tab
