#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates rNMP frequencies (3' position of aligned reads) and dNTPs located +/- 100 base pairs from rNMPs

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 5_Ribonucleotide-Frequencies.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of genome (sacCer2, nuclear, chrM, etc.)
	-r Reference genome assembly version (sacCer2, etc.)
	-d Local directory (/projects/home/agombolay3/data/repository/Ribose-seq-Project)"
}

#Use getopts function to create the command-line options ([-i], [-s], [-r], [-d], and [-h])
while getopts "i:s:r:d:h" opt; do
    case $opt in
        #Specify input as arrays to allow multiple input arguments
        i ) sample=($OPTARG) ;;
	#Specify input as variable to allow only one input argument
	s ) subset=$OPTARG ;;
	r ) reference=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #If user specifies [-h], print usage statement
        h ) usage ;;
    esac
done

#Exit program if user specifies [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Calculate nucleotide frequencies for each sample
for sample in ${sample[@]}; do

#############################################################################################################################
	#Input/Output
	
	#Location of input files
	referenceBED=$directory/ribose-seq/reference/$reference.bed
	referenceFasta1=$directory/ribose-seq/reference/$reference.fa
	#referenceFasta2=$directory/ribose-seq/reference/$reference.$subset.fa
	
	reads=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset/$sample.read-information.bed
	coordinates=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset/$sample.rNMP-coordinates.bed

	#Location of output directories
	output1=$directory/ribose-seq/results/Background-Frequencies
	
	output2=$directory/ribose-seq/results/$reference/$sample/Frequencies/rNMPs/$subset
	output3=$directory/ribose-seq/results/$reference/$sample/Frequencies/dNTPs/$subset/Results
	
	output4=$directory/ribose-seq/results/$reference/$sample/Frequencies/dNTPs/$subset/Columns/upstream
	output5=$directory/ribose-seq/results/$reference/$sample/Frequencies/dNTPs/$subset/Columns/downstream
	
	output6=$directory/ribose-seq/results/$reference/$sample/Frequencies/dNTPs/$subset/Columns/upstream/sequences
	output7=$directory/ribose-seq/results/$reference/$sample/Frequencies/dNTPs/$subset/Columns/downstream/sequences
	
	output8=$directory/ribose-seq/results/$reference/$sample/Frequencies/Datasets/$subset

	#Create directories if they do not already exist
	mkdir -p $output{1..8}

	#Remove any older versions of the output files
	rm -f $output{1..8}/*.txt
	
	#Location of output files
	background=$output1/Background-Frequencies.$reference.$subset.txt
	riboSequences=$output2/$sample.rNMP-Sequences.$reference.$subset.txt
	
	upstreamSequences=$output3/$sample.upstream-sequences.$reference.$subset.fa
	upstreamIntervals=$output3/$sample.upstream-intervals.$reference.$subset.bed
	
	downstreamSequences=$output3/$sample.downstream-sequences.$reference.$subset.fa
	downstreamIntervals=$output3/$sample.downstream-intervals.$reference.$subset.bed
	
	sequences1=$output6/$sample.upstream-sequences.$reference.$subset.txt
	sequences2=$output7/$sample.downstream-sequences.$reference.$subset.txt
	
	columns1=$output6/$sample.upstream-sequences.$reference.$subset.tab
	columns2=$output7/$sample.downstream-sequences.$reference.$subset.tab

	dataset=$output8/$sample.nucleotide-frequencies-dataset.$reference.$subset.txt
	zoomed=$output8/$sample.nucleotide-frequencies-zoomed.$reference.$subset.txt
		
#############################################################################################################################
	#STEP 1: Calculate background dNTP frequencies of reference genome

	#Index reference FASTA file
	#samtools faidx $referenceFasta1
	
	if [ $reference == "sacCer2" ] && [ $subset == "nuclear" ] ; then
		#Select only nuclear DNA and output to new file
		referenceFasta2=$(samtools faidx $referenceFasta1 | samtools faidx - 2micron chrI chrII \
		chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI)
	
	elif [ $reference == "sacCer2" ] && [ $subset == "chrM" ] ; then
		#Select only mitochondrial DNA and output to new file
		samtools faidx $referenceFasta1 chrM > $referenceFasta2
	fi
	
	#Calculate counts of each dNTP
	A_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'A' - | wc -l)
	C_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'C' - | wc -l)
	G_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'G' - | wc -l)
	T_Count0=$(grep -v '>' $referenceFasta2 | grep -o 'T' - | wc -l)
	
	#Calculate total number of dNTPs
	total0=$(($A_Count0+$C_Count0+$G_Count0+$T_Count0))

	#Calculate frequency of each dNTP
	A_Frequency0=$(echo "scale = 12; $A_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')
	C_Frequency0=$(echo "scale = 12; $C_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')
	G_Frequency0=$(echo "scale = 12; $G_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')
	T_Frequency0=$(echo "scale = 12; $T_Count0/$total0" | bc | awk '{printf "%.12f\n", $0}')

	#Save frequencies of dNTPs in the reference FASTA file (background frequencies) to TXT file
	echo -e "A\tC\tG\tU/T\n$A_Frequency0\t$C_Frequency0\t$G_Frequency0\t$U_Frequency0" > $background

#############################################################################################################################
	#STEP 2: Calculate rNMP Frequencies

	#Extract rNMP sequences
	if [ $subset == "nuclear" ]; then
		#Select only reads located in nuclear DNA and extract rNMP sequences
		grep -v 'chrM' $reads | awk '{print substr($0,length($0))}' - > $riboSequences
	elif [ $subset == "chrM" ]; then
		#Select only reads located in mitochondrial DNA and extract rNMP sequences
		grep 'chrM' $reads | awk '{print substr($0,length($0))}' - > $riboSequences
	else
		#Select all reads located in genomic DNA and extract rNMP sequences
		cat $reads | awk '{print substr($0,length($0))}' - > $riboSequences
	fi
	
	#Calculate counts of rNMPs
	A_Count1=$(awk '$1 == "A"' $riboSequences | wc -l)
	C_Count1=$(awk '$1 == "C"' $riboSequences | wc -l)
	G_Count1=$(awk '$1 == "G"' $riboSequences | wc -l)
	U_Count1=$(awk '$1 == "T"' $riboSequences | wc -l)
	
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
	#STEP 3: Obtain coordinates and sequences of +/- 100 downstream/upstream dNTPs from rNMPs

	#Obtain coordinates of upstream/downstream sequences based on rNMP coordinates
	bedtools flank -i $coordinates -s -g $referenceBED -l 100 -r 0 > $upstreamIntervals
	bedtools flank -i $coordinates -s -g $referenceBED -l 0 -r 100 > $downstreamIntervals

	#Obtain sequences of upstream/downstream coordinates
	bedtools getfasta -s -fi $referenceFasta2 -bed $upstreamIntervals -fo $upstreamSequences
	bedtools getfasta -s -fi $referenceFasta2 -bed $downstreamIntervals -fo $downstreamSequences

#############################################################################################################################
	#STEP 4: Tabulate sequences of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Extract sequences from FASTA files
	#Reverse order of upstream nucleotides
	grep -v '>' $upstreamSequences | rev > $sequences1
	grep -v '>' $downstreamSequences > $sequences2
				
	#Insert tabs between each nucleotide
	cat $sequences1 | sed 's/.../& /2g;s/./& /g' > $columns1
	cat $sequences2 | sed 's/.../& /2g;s/./& /g' > $columns2

	for i in {1..100}; do
		#Location of output files
		lists1=$output4/$sample.column.$i.upstream.$reference.$subset.txt
		lists2=$output5/$sample.column.$i.downstream.$reference.$subset.txt
		#Save lists of dNTPs at each +/- 100 bp downstream/upstream position
		awk -v field=$i '{ print $field }' $columns1 > $lists1
		awk -v field=$i '{ print $field }' $columns2 > $lists2
	done

#############################################################################################################################
	#STEP 5: Calculate frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Calculate frequencies at each position
	for file in `ls -v $output4/$sample*.txt`; do

		#Calculate count of each dNTP
		A_Count2=$(grep -o 'A' $file | wc -l)
		C_Count2=$(grep -o 'C' $file | wc -l)
		G_Count2=$(grep -o 'G' $file | wc -l)
		T_Count2=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNTPs
		total2=$(($A_Count2+$C_Count2+$G_Count2+$T_Count2))

		#Calculate normalized frequencies of dNTPs
		A_Frequency2=$(echo "scale = 12; ($A_Count2/$total2)/$A_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		C_Frequency2=$(echo "scale = 12; ($C_Count2/$total2)/$C_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		G_Frequency2=$(echo "scale = 12; ($G_Count2/$total2)/$G_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		T_Frequency2=$(echo "scale = 12; ($T_Count2/$total2)/$T_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNTPs frequencies to TXT file
		echo $A_Frequency2 >> A_frequency2.txt; echo $C_Frequency2 >> C_frequency2.txt
		echo $G_Frequency2 >> G_frequency2.txt; echo $T_Frequency2 >> T_frequency2.txt
			
		#Save upstream dNTP frequencies together and reverse order of frequencies so ordered from -100 --> -1
		upstreamFrequencies=$(paste A_frequency2.txt C_frequency2.txt G_frequency2.txt T_frequency2.txt | tac -)
		
	done

	#Calculate frequencies at each position
	for file in `ls -v $output5/$sample*.txt`; do

		#Calculate count of each dNTP
		A_Count3=$(grep -v '>' $file | grep -o 'A' - | wc -l)
		C_Count3=$(grep -v '>' $file | grep -o 'C' - | wc -l)
		G_Count3=$(grep -v '>' $file | grep -o 'G' - | wc -l)
		T_Count3=$(grep -v '>' $file | grep -o 'T' - | wc -l)

		#Calculate total number of dNTPs
		total3=$(($A_Count3+$C_Count3+$G_Count3+$T_Count3))
	
		#Calculate normalized frequencies of dNTPs
		A_Frequency3=$(echo "scale = 12; ($A_Count3/$total3)/$A_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		C_Frequency3=$(echo "scale = 12; ($C_Count3/$total3)/$C_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		G_Frequency3=$(echo "scale = 12; ($G_Count3/$total3)/$G_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		T_Frequency3=$(echo "scale = 12; ($T_Count3/$total3)/$T_Frequency0" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized dNTPs frequencies to TXT file
		echo $A_Frequency3 >> A_frequency3.txt; echo $C_Frequency3 >> C_frequency3.txt
		echo $G_Frequency3 >> G_frequency3.txt; echo $T_Frequency3 >> T_frequency3.txt
		
		#Save upstream dNTP frequencies together (frequencies will be ordered from +1 --> +100)
		downstreamFrequencies=$(paste A_frequency3.txt C_frequency3.txt G_frequency3.txt T_frequency3.txt)
	
	done

	#Remove intermediate files
	rm A_frequency{2..3}.txt C_frequency{2..3}.txt G_frequency{2..3}.txt T_frequency{2..3}.txt
	
#############################################################################################################################
	#STEP 6: Create dataset file containing nucleotide frequencies for plotting

	#Combine rNMP frequencies and upstream and downstream dNTP frequencies in appropriate order
	data1=$(cat <(echo "$upstreamFrequencies") <(echo "$riboFrequencies") <(echo "$downstreamFrequencies"))
	
	#Add nucleotide positions (-100 --> +100) and add nucleotides to header line (A, C, G, and U/T)
	echo -e "\tA\tC\tG\tU/T" > $dataset && paste <(echo "$(seq -100 1 100)") <(cat <(echo "$data1")) >> $dataset

	#Smaller dataset (-15 nt to +15 nt)
	data2=$(head -117 $dataset | tail -31)
	echo -e "\tA\tC\tG\tU/T" > $zoomed && cat <(echo "$data2") >> $zoomed

	#Let the user know the program is has finished running
	echo "Calculation of nucleotide frequencies for $sample is complete"

done
