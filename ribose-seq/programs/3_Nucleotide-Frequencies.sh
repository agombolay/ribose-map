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
	
	#Location of "Reference" directory
	directory0=$directory/ribose-seq/reference/

	#Location of "Nucleotide-Frequencies" directory
	directory2=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies

##########################################################################################################################################
	#STEP 1: Covert BAM alignment file to FASTA format

	#Location of input file
	bam=$directory1/$sample.bam
	echo $bam

	#Location of output directory
	output1=$directory2/rNMPs/$subset

	#Create directory if it does not already exist
	mkdir -p $output1

	#Remove previously created files so new files are created
	rm -f $output1/*.txt
	
	#Location of output files	
	fastq=$output1/$sample.aligned-reads.fq
	fasta=$output1/$sample.aligned-reads.fa

	#Convert BAM file to FASTQ
	samtools bam2fq $bam > $fastq

	#Convert FASTQ file to FASTA
	seqtk seq -A $fastq > $fasta
	
	#Output only sequences in FASTA file (exclude all header lines)
	grep -v '>' $fasta > temporary && mv temporary $fasta
	
##########################################################################################################################################
	#STEP 2: Obtain rNMP coordinates from aligned reads

	#Location of output files
	bed=$output1/$sample.aligned-reads.bed
	coverage=$output1/$sample.rNMP-coverage.0-based.txt
	
	readCoordinates=$output1/$sample.read-coordinates.bed
	readInformation=$output1/$sample.read-information.bed
	
	riboCoordinates1=$output1/$sample.rNMP-coordinates.0-based.genome.bed
	riboCoordinates2=$output1/$sample.rNMP-coordinates.0-based.$subset.bed
	
	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > $bed
	
	#Obtain genome coverage of aligned reads
	bedtools genomecov -3 -bg -ibam $bam > $coverage
	
	#Extract aligned read coordinates, sequences, and strands from BED and SAM files
	paste $bed $fasta | awk -v "OFS=\t" '{print $1, $2, $3, $4, $6, $7}' > $readInformation
	
	#Determine rNMP coordinates from reads aligned to positive strand
	awk -v "OFS=\t" '$5 == "+" {print $1, ($3 - 1), $3, " ", " ", $5}' $readInformation > positive-reads.txt
	
	#Determine rNMP coordinates from reads aligned to negative strand
	awk -v "OFS=\t" '$5 == "-" {print $1, $2, ($2 + 1), " ", " ", $5}' $readInformation > negative-reads.txt
	cat positive-reads.txt negative-reads.txt > $riboCoordinates1
	
	#Select only rNMP coordinates located in nuclear DNA
	if [ $subset == "nuclear" ]; then
		grep -v 'chrM' $riboCoordinates1 > $riboCoordinates2
	#Select only rNMP coordinates located in mitochondrial DNA
	elif [ $subset == "chrM" ]; then
		grep 'chrM' $riboCoordinates1 > $riboCoordinates2
	#Select all rNMP coordinates located in genomic DNA
	else
		cat $riboCoordinates1 > $riboCoordinates2
	fi
	
	#Remove intermediate files
	rm positive-reads.txt negative-reads.txt
##########################################################################################################################################
	#STEP 3: Calculate background dNTP frequencies of reference genome

	#Location of input file
	referenceFasta1=$directory0/$subset.fa

	#Location of output directory
	output2=$directory/ribose-seq/results/Background-dNTP-Frequencies

	#Create directory if it does not already exist
	mkdir -p $output2

	#Location of output file
	background=$output2/$reference.$subset.Background-dNTP-Frequencies.txt

	#Remove previously created files so new files are created
	rm -f $output2/*.txt
	
	#Calculate counts of each dNTP
	A_backgroundCount=$(grep -v '>' $referenceFasta1 | grep -o 'A' - | wc -l)
	C_backgroundCount=$(grep -v '>' $referenceFasta1 | grep -o 'C' - | wc -l)
	G_backgroundCount=$(grep -v '>' $referenceFasta1 | grep -o 'G' - | wc -l)
	T_backgroundCount=$(grep -v '>' $referenceFasta1 | grep -o 'T' - | wc -l)
	
	#Calculate total number of dNTPs
	total_backgroundCount=$(($A_backgroundCount+$C_backgroundCount+$G_backgroundCount+$T_backgroundCount))

	#Calculate frequency of each dNTP
	A_backgroundFrequency=$(echo "scale = 12; $A_backgroundCount/$total_backgroundCount" | bc | awk '{printf "%.12f\n", $0}')
	C_backgroundFrequency=$(echo "scale = 12; $C_backgroundCount/$total_backgroundCount" | bc | awk '{printf "%.12f\n", $0}')
	G_backgroundFrequency=$(echo "scale = 12; $G_backgroundCount/$total_backgroundCount" | bc | awk '{printf "%.12f\n", $0}')
	T_backgroundFrequency=$(echo "scale = 12; $T_backgroundCount/$total_backgroundCount" | bc | awk '{printf "%.12f\n", $0}')

	#Save frequencies of dNTPs to TXT file
	echo "A Background Frequency: $A_backgroundFrequency" >> $background
	echo "C Background Frequency: $C_backgroundFrequency" >> $background
	echo "G Background Frequency: $G_backgroundFrequency" >> $background
	echo "T Background Frequency: $T_backgroundFrequency" >> $background

##########################################################################################################################################
	#STEP 4: Calculate rNMP Frequencies

	#Location of output files
	riboSequences=$output1/$sample.rNMP-Sequences.$reference.$subset.txt
	riboFrequencies=$output1/$sample.rNMP-frequencies.$reference.$subset.txt	
	
	#Select only reads located in nuclear DNA
	if [ $subset == "nuclear" ]; then
		grep -v 'chrM' $readInformation > temporary
	#Select only reads located in mitochondrial DNA
	elif [ $subset == "chrM" ]; then
		grep 'chrM' $readInformation > temporary
	#Select all reads located in genomic DNA
	else
		cat $readInformation > temporary
	fi
	
	#Extract rNMP sequences from 3' end of aligned reads
	awk '{print substr($0,length($0))}' temporary > $riboSequences
	
	#Calculate counts of rNMPs
	A_riboCount=$(awk '$1 == "A"' $riboSequences | wc -l)
	C_riboCount=$(awk '$1 == "C"' $riboSequences | wc -l)
	G_riboCount=$(awk '$1 == "G"' $riboSequences | wc -l)
	U_riboCount=$(awk '$1 == "T"' $riboSequences | wc -l)
	
	#Calculate total number of rNMPs
	total_riboCount=$(($A_riboCount+$C_riboCount+$G_riboCount+$U_riboCount))

	#Calculate raw frequency of each rNMP
	A_rawRiboFrequency=$(echo "scale = 12; $A_riboCount/$total_riboCount" | bc | awk '{printf "%.12f\n", $0}')
	C_rawRiboFrequency=$(echo "scale = 12; $C_riboCount/$total_riboCount" | bc | awk '{printf "%.12f\n", $0}')
	G_rawRiboFrequency=$(echo "scale = 12; $G_riboCount/$total_riboCount" | bc | awk '{printf "%.12f\n", $0}')
	U_rawRiboFrequency=$(echo "scale = 12; $U_riboCount/$total_riboCount" | bc | awk '{printf "%.12f\n", $0}')

	#Calculate normalized frequency of each rNMP
	A_riboFrequency=$(echo "scale = 12; $A_rawRiboFrequency/$A_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
	C_riboFrequency=$(echo "scale = 12; $C_rawRiboFrequency/$C_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
	G_riboFrequency=$(echo "scale = 12; $G_rawRiboFrequency/$G_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
	U_riboFrequency=$(echo "scale = 12; $U_rawRiboFrequency/$T_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')

	#Save normalized frequencies of rNMPs to TXT file
	echo -e "$A_riboFrequency\t$C_riboFrequency\t$G_riboFrequency\t$U_riboFrequency" > $riboFrequencies

	rm temporary
##########################################################################################################################################
	#STEP 5: Obtain coordinates and sequences of +/- 100 downstream/upstream dNTPs from rNMPs

	#Location of input files
	referenceBED=$directory0/$reference.bed
	referenceFasta2=$directory0/$reference.fa

	#Location of output directory
	output3=$directory2/dNTPs/$subset

	#Create directory if it does not already exist
    	mkdir -p $output3

	#Location of output files
	upstreamSequences=$output3/$sample.upstream-sequences.$reference.$subset.fa
	upstreamIntervals=$output3/$sample.upstream-intervals.$reference.$subset.bed
	downstreamSequences=$output3/$sample.downstream-sequences.$reference.$subset.fa
	downstreamIntervals=$output3/$sample.downstream-intervals.$reference.$subset.bed
	
	#Obtain coordinates of upstream/downstream sequences based on rNMP coordinates
	bedtools flank -i $riboCoordinates2 -s -g $referenceBED -l 100 -r 0 > $upstreamIntervals
	bedtools flank -i $riboCoordinates2 -s -g $referenceBED -l 0 -r 100 > $downstreamIntervals

	#Obtain sequences of upstream/downstream coordinates
	bedtools getfasta -s -fi $referenceFasta2 -bed $upstreamIntervals -fo $upstreamSequences
	bedtools getfasta -s -fi $referenceFasta2 -bed $downstreamIntervals -fo $downstreamSequences

##########################################################################################################################################
	#STEP 6: Tabulate sequences of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Location of output directory
	output4=$directory2/dNTPs/$subset/Columns/upstream
	output5=$directory2/dNTPs/$subset/Columns/downstream
	output6=$directory2/dNTPs/$subset/Columns/upstream/sequences
	output7=$directory2/dNTPs/$subset/Columns/downstream/sequences
				
	#Create directories if they do not already exist
    	mkdir -p $output4 $output5 $output6 $output7
	
	#Remove previously created files so new files are created
	rm -f $output4/*.txt $output5/*.txt $output6/*.txt $output7/*.txt
			
	sequences1=$output6/$sample.upstream-sequences.$reference.$subset.txt
	sequences2=$output7/$sample.downstream-sequences.$reference.$subset.txt
	sequences1Reversed=$output6/$sample.upstream-sequences.$reference.$subset.reversed.txt
	
	grep -v '>' $upstreamSequences > $sequences1
	grep -v '>' $downstreamSequences > $sequences2
	
	cat $sequences1|rev > $sequences1Reversed
	
	columns1=$output6/$sample.upstream-sequences.$reference.$subset.tab
	columns2=$output7/$sample.downstream-sequences.$reference.$subset.tab
				
	#Insert tabs between each nucleotide
	cat $sequences1Reversed | sed 's/.../& /2g;s/./& /g' > $columns1
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
	#STEP 7: Calculate frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Location of output directory
	output8=$directory2/dNTPs/$subset/Raw-Data

	#Create directory if it does not already exist
	mkdir -p $output8

	#Remove previously created files so new files are created
	rm -f $output8/*.txt
	
	#Location of output files (indivdiual base frequencies)
	A_upstreamFrequencies=$output8/A_dNTP-frequencies.$reference.$subset.upstream.txt
	C_upstreamFrequencies=$output8/C_dNTP-frequencies.$reference.$subset.upstream.txt
	G_upstreamFrequencies=$output8/G_dNTP-frequencies.$reference.$subset.upstream.txt
	T_upstreamFrequencies=$output8/T_dNTP-frequencies.$reference.$subset.upstream.txt

	#Location of output file (combined base frequencies)
	upstreamFrequencies=$output8/$sample.dNTP-frequencies.$reference.$subset.upstream.txt
	
	#Calculate dNTP frequencies for each 100 upstream position
	for file in `ls -v $directory2/dNTPs/$subset/Columns/upstream/$sample*.txt`; do

		#Calculate count of each dNTP
		A_upstreamCount=$(grep -o 'A' $file | wc -l)
		C_upstreamCount=$(grep -o 'C' $file | wc -l)
		G_upstreamCount=$(grep -o 'G' $file | wc -l)
		T_upstreamCount=$(grep -o 'T' $file | wc -l)

		#Calculate total number of dNTPs
		total_upstreamCount=$(($A_upstreamCount+$C_upstreamCount+$G_upstreamCount+$T_upstreamCount))
	
		#Calculate raw frequencies of dNTPs
		A_rawUpstreamFrequency=$(echo "scale = 12; $A_upstreamCount/$total_upstreamCount" | bc | awk '{printf "%.12f\n", $0}')
		C_rawUpstreamFrequency=$(echo "scale = 12; $C_upstreamCount/$total_upstreamCount" | bc | awk '{printf "%.12f\n", $0}')
		G_rawUpstreamFrequency=$(echo "scale = 12; $G_upstreamCount/$total_upstreamCount" | bc | awk '{printf "%.12f\n", $0}')
		T_rawUpstreamFrequency=$(echo "scale = 12; $T_upstreamCount/$total_upstreamCount" | bc | awk '{printf "%.12f\n", $0}')

		#Calculate normalized frequencies of dNTPs
		A_upstreamFrequency=$(echo "scale = 12; $A_rawUpstreamFrequency/$A_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		C_upstreamFrequency=$(echo "scale = 12; $C_rawUpstreamFrequency/$C_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		G_upstreamFrequency=$(echo "scale = 12; $G_rawUpstreamFrequency/$G_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		T_upstreamFrequency=$(echo "scale = 12; $T_rawUpstreamFrequency/$T_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized frequencies of dNTPs to TXT file
		echo $A_upstreamFrequency >> $A_upstreamFrequencies
		echo $C_upstreamFrequency >> $C_upstreamFrequencies
		echo $G_upstreamFrequency >> $G_upstreamFrequencies
		echo $T_upstreamFrequency >> $T_upstreamFrequencies

		#Save frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs to one TXT file
		paste $A_upstreamFrequencies $C_upstreamFrequencies $G_upstreamFrequencies $T_upstreamFrequencies > $upstreamFrequencies
	
		#Reverse order of frequencies so ordered from -100 --> -1
		tac $upstreamFrequencies > temporary && mv temporary $upstreamFrequencies
		
	done
	
	#Location of output files (indivdiual base frequencies)
	A_downstreamFrequencies=$output8/A_dNTP-frequencies.$reference.$subset.downstream.txt
	C_downstreamFrequencies=$output8/C_dNTP-frequencies.$reference.$subset.downstream.txt
	G_downstreamFrequencies=$output8/G_dNTP-frequencies.$reference.$subset.downstream.txt
	T_downstreamFrequencies=$output8/T_dNTP-frequencies.$reference.$subset.downstream.txt

	#Location of output file (combined base frequencies)
	downstreamFrequencies=$output8/$sample.dNTP-frequencies.$reference.$subset.downstream.txt
	
	#Calculate dNTP frequencies for each 100 downstream position
	for file in `ls -v $directory2/dNTPs/$subset/Columns/downstream/$sample*.txt`; do

		#Calculate count of each dNTP
		A_downstreamCount=$(grep -v '>' $file | grep -o 'A' - | wc -l)
		C_downstreamCount=$(grep -v '>' $file | grep -o 'C' - | wc -l)
		G_downstreamCount=$(grep -v '>' $file | grep -o 'G' - | wc -l)
		T_downstreamCount=$(grep -v '>' $file | grep -o 'T' - | wc -l)

		#Calculate total number of dNTPs
		total_downstreamCount=$(($A_downstreamCount+$C_downstreamCount+$G_downstreamCount+$T_downstreamCount))
	
		#Calculate raw frequencies of dNTPs
		A_rawDownstreamFrequency=$(echo "scale = 12; $A_downstreamCount/$total_downstreamCount" | bc | awk '{printf "%.12f\n", $0}')
		C_rawDownstreamFrequency=$(echo "scale = 12; $C_downstreamCount/$total_downstreamCount" | bc | awk '{printf "%.12f\n", $0}')
		G_rawDownstreamFrequency=$(echo "scale = 12; $G_downstreamCount/$total_downstreamCount" | bc | awk '{printf "%.12f\n", $0}')
		T_rawDownstreamFrequency=$(echo "scale = 12; $T_downstreamCount/$total_downstreamCount" | bc | awk '{printf "%.12f\n", $0}')

		#Calculate normalized frequencies of dNTPs
		A_downstreamFrequency=$(echo "scale = 12; $A_rawDownstreamFrequency/$A_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		C_downstreamFrequency=$(echo "scale = 12; $C_rawDownstreamFrequency/$C_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		G_downstreamFrequency=$(echo "scale = 12; $G_rawDownstreamFrequency/$G_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		T_downstreamFrequency=$(echo "scale = 12; $T_rawDownstreamFrequency/$T_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		
		#Save normalized frequencies of dNTPs to TXT file
		echo $A_downstreamFrequency >> $A_downstreamFrequencies
		echo $C_downstreamFrequency >> $C_downstreamFrequencies
		echo $G_downstreamFrequency >> $G_downstreamFrequencies
		echo $T_downstreamFrequency >> $T_downstreamFrequencies

		#Save frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs to one TXT file
		paste $A_downstreamFrequencies $C_downstreamFrequencies $G_downstreamFrequencies $T_downstreamFrequencies > $downstreamFrequencies
	done

##########################################################################################################################################
	#STEP 8: Create dataset file containing nucleotide frequencies needed for plotting

	#Location of output directory
	output9=$directory2/Datasets/$subset

	#Create directory if it does not already exist
    	mkdir -p $output9

	#Location of output file
	dataset=$output9/$sample.nucleotide-frequencies-dataset.$reference.$subset.txt
	zoomed=$output9/$sample.nucleotide-frequencies-zoomed.$reference.$subset.txt
	
	#Remove previously created files so new ones are created
	rm -f $output9/*.txt
	
	#Print values -100 to 100
	seq -100 1 100 > temporary1
	
	#Save files containing rNMP and upstream/downstream dNTP frequencies to one file
	cat $upstreamFrequencies $riboFrequencies $downstreamFrequencies >> temporary2

	#Save files to one combined TXT file
	paste temporary1 temporary2 > temporary3

	#Add header line containing nucleotides to beginning of file 
	echo -e "\tA\tC\tG\tU/T" > $dataset; cat temporary3 >> $dataset;
	
	#Smaller dataset
	head -117 $dataset | tail -31 > temporary4
	echo -e "\tA\tC\tG\tU/T" > $zoomed; cat temporary4 >> $zoomed;

	#Remove temporary files
	rm -f temporary1 temporary2 temporary3 temporary4

	echo "Calculation of nucleotide frequencies for $sample is complete"

done
