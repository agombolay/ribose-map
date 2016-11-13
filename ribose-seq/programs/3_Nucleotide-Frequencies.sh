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

##########################################################################################################################################
	#Input/Output
	
	#Location of input files
	referenceBED=$directory/ribose-seq/reference/$reference.bed
	referenceFasta=$directory/ribose-seq/reference/$reference.fa
	
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
	riboFrequencies=$output2/$sample.rNMP-frequencies.$reference.$subset.txt
	
	upstreamSequences=$output3/$sample.upstream-sequences.$reference.$subset.fa
	upstreamIntervals=$output3/$sample.upstream-intervals.$reference.$subset.bed
	
	downstreamSequences=$output3/$sample.downstream-sequences.$reference.$subset.fa
	downstreamIntervals=$output3/$sample.downstream-intervals.$reference.$subset.bed
	
	sequences1=$output6/$sample.upstream-sequences.$reference.$subset.txt
	sequences2=$output7/$sample.downstream-sequences.$reference.$subset.txt
	reversed=$output6/$sample.upstream-sequences.$reference.$subset.reversed.txt
	
	columns1=$output6/$sample.upstream-sequences.$reference.$subset.tab
	columns2=$output7/$sample.downstream-sequences.$reference.$subset.tab
		
	upstreamFrequencies=$output3/$sample.dNTP-frequencies.$reference.$subset.upstream.txt
	downstreamFrequencies=$output3/$sample.dNTP-frequencies.$reference.$subset.downstream.txt
	
	dataset=$output8/$sample.nucleotide-frequencies-dataset.$reference.$subset.txt
	zoomed=$output8/$sample.nucleotide-frequencies-zoomed.$reference.$subset.txt
		
##########################################################################################################################################
	#STEP 1: Calculate background dNTP frequencies of reference genome

	#Index reference FASTA file
	samtools faidx $referenceFasta
	
	if [ $reference == "sacCer2" ] && [ $subset == "nuclear" ] ; then
		#Select only nuclear DNA and output to new file
		samtools faidx $referenceFasta 2micron chrI chrII chrIII chrIV chrV chrVI chrVII \
		chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI > sacCer2.nuclear.fa
	
	elif [ $reference == "sacCer2" ] && [ $subset == "chrM" ] ; then
		#Select only mitochondrial DNA and output to new file
		samtools faidx $referenceFasta chrM > sacCer2.chrM.fa
	fi
	
	#Calculate counts of each dNTP
	A_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'A' - | wc -l)
	C_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'C' - | wc -l)
	G_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'G' - | wc -l)
	T_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'T' - | wc -l)
	
	#Calculate total number of dNTPs
	total1=$(($A_backgroundCount+$C_backgroundCount+$G_backgroundCount+$T_backgroundCount))

	#Calculate frequency of each dNTP
	A_backgroundFrequency=$(echo "scale = 12; $A_backgroundCount/$total1" | bc | awk '{printf "%.12f\n", $0}')
	C_backgroundFrequency=$(echo "scale = 12; $C_backgroundCount/$total1" | bc | awk '{printf "%.12f\n", $0}')
	G_backgroundFrequency=$(echo "scale = 12; $G_backgroundCount/$total1" | bc | awk '{printf "%.12f\n", $0}')
	T_backgroundFrequency=$(echo "scale = 12; $T_backgroundCount/$total1" | bc | awk '{printf "%.12f\n", $0}')

	#Save frequencies of background dNTPs in the reference FASTA file to TXT file
	echo -e "A\tC\tG\tU/T\n$A_backgroundFrequency\t$C_backgroundFrequency\t$G_backgroundFrequency\t$U_backgroundFrequency" \
	> $background

##########################################################################################################################################
	#STEP 2: Calculate rNMP Frequencies

	#Select only reads located in nuclear DNA
	if [ $subset == "nuclear" ]; then
		grep -v 'chrM' $reads > temporary
	#Select only reads located in mitochondrial DNA
	elif [ $subset == "chrM" ]; then
		grep 'chrM' $reads > temporary
	#Select all reads located in genomic DNA
	else
		cat $reads > temporary
	fi
	
	#Extract rNMP sequences from 3' end of aligned reads
	awk '{print substr($0,length($0))}' temporary > $riboSequences
	
	#Calculate counts of rNMPs
	A_riboCount=$(awk '$1 == "A"' $riboSequences | wc -l)
	C_riboCount=$(awk '$1 == "C"' $riboSequences | wc -l)
	G_riboCount=$(awk '$1 == "G"' $riboSequences | wc -l)
	U_riboCount=$(awk '$1 == "T"' $riboSequences | wc -l)
	
	#Calculate total number of rNMPs
	total2=$(($A_riboCount+$C_riboCount+$G_riboCount+$U_riboCount))

	#Calculate normalized frequency of each rNMP
	A_riboFrequency=$(echo "scale = 12; ($A_riboCount/$total2)/$A_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
	C_riboFrequency=$(echo "scale = 12; ($C_riboCount/$total2)/$C_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
	G_riboFrequency=$(echo "scale = 12; ($G_riboCount/$total2)/$G_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
	U_riboFrequency=$(echo "scale = 12; ($U_riboCount/$total2)/$T_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')

	#Save normalized frequencies of rNMPs to TXT file
	echo -e "$A_riboFrequency\t$C_riboFrequency\t$G_riboFrequency\t$U_riboFrequency" > $riboFrequencies

	#Remove temporary file
	rm temporary
	
##########################################################################################################################################
	#STEP 3: Obtain coordinates and sequences of +/- 100 downstream/upstream dNTPs from rNMPs

	#Obtain coordinates of upstream/downstream sequences based on rNMP coordinates
	bedtools flank -i $riboCoordinates2 -s -g $referenceBED -l 100 -r 0 > $upstreamIntervals
	bedtools flank -i $riboCoordinates2 -s -g $referenceBED -l 0 -r 100 > $downstreamIntervals

	#Obtain sequences of upstream/downstream coordinates
	bedtools getfasta -s -fi $referenceFasta -bed $upstreamIntervals -fo $upstreamSequences
	bedtools getfasta -s -fi $referenceFasta -bed $downstreamIntervals -fo $downstreamSequences

##########################################################################################################################################
	#STEP 4: Tabulate sequences of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Extract sequences from FASTA files
	grep -v '>' $upstreamSequences > $sequences1
	grep -v '>' $downstreamSequences > $sequences2
	
	#Reverse order of upstream nucleotides
	cat $sequences1|rev > $reversed
				
	#Insert tabs between each nucleotide
	cat $reversed | sed 's/.../& /2g;s/./& /g' > $columns1
	cat $sequences2 | sed 's/.../& /2g;s/./& /g' > $columns2

	for i in {1..100}; do
		#Location of output files
		lists1=$output4/$sample.column.$i.upstream.$reference.$subset.txt
		lists2=$output5/$sample.column.$i.downstream.$reference.$subset.txt
		#Save lists of dNTPs at each +/- 100 bp downstream/upstream position
		awk -v field=$i '{ print $field }' $columns1 > $lists1
		awk -v field=$i '{ print $field }' $columns2 > $lists2
	done

##########################################################################################################################################
	#STEP 5: Calculate frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Calculate frequencies at each position
	for file in `ls -v $output4/$sample*.txt`; do

		#Calculate count of each dNTP
		A_upstreamCount=$(grep -o 'A' $file | wc -l)
		C_upstreamCount=$(grep -o 'C' $file | wc -l)
		G_upstreamCount=$(grep -o 'G' $file | wc -l)
		T_upstreamCount=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNTPs
		total3=$(($A_upstreamCount+$C_upstreamCount+$G_upstreamCount+$T_upstreamCount))

		#Calculate normalized frequencies of dNTPs
		A_upstreamFrequency=$(echo "scale = 12; ($A_upstreamCount/$total3)/$A_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		C_upstreamFrequency=$(echo "scale = 12; ($C_upstreamCount/$total3)/$C_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		G_upstreamFrequency=$(echo "scale = 12; ($G_upstreamCount/$total3)/$G_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		T_upstreamFrequency=$(echo "scale = 12; ($T_upstreamCount/$total3)/$T_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized frequencies of dNTPs to TXT file
		echo $A_upstreamFrequency >> A_frequencies1.txt; echo $C_upstreamFrequency >> C_frequencies1.txt
		echo $G_upstreamFrequency >> G_frequencies1.txt; echo $T_upstreamFrequency >> T_frequencies1.txt
			
		#Save frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs to one TXT file
		paste A_frequencies1.txt C_frequencies1.txt G_frequencies1.txt T_frequencies1.txt > $upstreamFrequencies
	
		#Reverse order of nucleotide frequencies so ordered from -100 --> -1
		tac $upstreamFrequencies > temporary && mv temporary $upstreamFrequencies
		
	done
	
	#Calculate frequencies at each position
	for file in `ls -v $output5/$sample*.txt`; do

		#Calculate count of each dNTP
		A_downstreamCount=$(grep -v '>' $file | grep -o 'A' - | wc -l)
		C_downstreamCount=$(grep -v '>' $file | grep -o 'C' - | wc -l)
		G_downstreamCount=$(grep -v '>' $file | grep -o 'G' - | wc -l)
		T_downstreamCount=$(grep -v '>' $file | grep -o 'T' - | wc -l)

		#Calculate total number of dNTPs
		total4=$(($A_downstreamCount+$C_downstreamCount+$G_downstreamCount+$T_downstreamCount))
	
		#Calculate normalized frequencies of dNTPs
		A_downstreamFrequency=$(echo "scale = 12; ($A_downstreamCount/$total4)/$A_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		C_downstreamFrequency=$(echo "scale = 12; ($C_downstreamCount/$total4)/$C_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		G_downstreamFrequency=$(echo "scale = 12; ($G_downstreamCount/$total4)/$G_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		T_downstreamFrequency=$(echo "scale = 12; ($T_downstreamCount/$total4)/$T_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized frequencies of dNTPs to TXT file
		echo $A_downstreamFrequency >> A_frequencies2.txt; echo $C_downstreamFrequency >> C_frequencies2.txt
		echo $G_downstreamFrequency >> G_frequencies2.txt; echo $T_downstreamFrequency >> T_frequencies2.txt
		
		#Save frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs to one TXT file
		paste A_frequencies2.txt C_frequencies2.txt G_frequencies2.txt T_frequencies2.txt > $downstreamFrequencies
	
	done

	#Remove intermediate files
	rm A_frequencies{1..2}.txt C_frequencies{1..2}.txt G_frequencies{1..2}.txt T_frequencies{1..2}.txt
	
##########################################################################################################################################
	#STEP 6: Create dataset file containing nucleotide frequencies needed for plotting

	#Print values -100 to 100
	seq -100 1 100 > temporary1
	
	#Save files containing rNMP and upstream/downstream dNTP frequencies to file
	cat $upstreamFrequencies $riboFrequencies $downstreamFrequencies >> temporary2

	#Save files to one combined TXT file
	paste temporary1 temporary2 > temporary3

	#Add header line containing nucleotides to beginning of file 
	echo -e "\tA\tC\tG\tU/T" > $dataset; cat temporary3 >> $dataset;
	
	#Smaller dataset (-15 nt to +15 nt)
	head -117 $dataset | tail -31 > temporary4
	echo -e "\tA\tC\tG\tU/T" > $zoomed; cat temporary4 >> $zoomed;

	#Remove temporary files
	rm -f temporary1 temporary2 temporary3 temporary4

	#Let the user know the program is has finished running
	echo "Calculation of nucleotide frequencies for $sample is complete"

done
