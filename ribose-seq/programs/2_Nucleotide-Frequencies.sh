#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates rNMP frequencies (3' position of aligned reads) and dNTPs located +/- 100 base pairs from rNMPs

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 5_Ribonucleotide-Frequencies.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of genome (sacCer2, nuclear, mitochondria)
	-r Reference genome assembly version (sacCer2, etc.)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
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

	#Location of "Alignment" directory
	directory1=$directory/ribose-seq/results/$reference/$sample/Alignment

	#Location of "Nucleotide-Frequencies" directory
	directory2=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies

##########################################################################################################################################
	#STEP 1: Covert BAM alignment file to FASTA format

	#Location of input file
	bam=$directory1/$sample.bam

	#Location of output directory
	output1=$directory2/rNMPs/$subset

	#Create directory if it does not already exist
	if [[ ! -d $output1 ]]; then
		mkdir -p $output1
	fi

	#Location of output files	
	fastq=$output1/$sample.aligned-reads.fastq
	fasta=$output1/$sample.aligned-reads.fasta

	#Convert BAM file to FASTQ
	samtools bam2fq $bam > $fastq

	#Convert FASTQ file to FASTA
	seqtk seq -A $fastq > $fasta

##########################################################################################################################################
	#STEP 2: Obtain rNMP coordinates from aligned reads

	#Location of output files
	bed=$output1/$sample.aligned-reads.bed
	sam=$output1/$sample.aligned-reads.sam
	readCoordinates=$output1/$sample.read-coordinates.bed
	positiveCoordinates0=$output1/$sample.rNMP-coordinates.positive.0-based.txt
	negativeCoordinates0=$output1/$sample.rNMP-coordinates.negative.0-based.txt
	positiveCoordinates1=$output1/$sample.rNMP-coordinates.positive.1-based.txt
	negativeCoordinates1=$output1/$sample.rNMP-coordinates.negative.1-based.txt
	
	#0-BASED COORDINATES of READS:
	#Covert BAM file to BED format
	bedtools bamtobed -i $bam > $bed

	#Convert BAM file to SAM format
	samtools view $bam > $sam

	#Extract aligned read coordinates, sequences, and strands from BED and SAM files
	paste $bed $sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' > $readCoordinates

	#0-BASED COORDINATES OF rNMPs:
	#Obtain coordinates of rNMPs (3’ end of aligned read):
	bedtools genomecov -3 -strand + -bg -ibam $bam > $positiveCoordinates0
	bedtools genomecov -3 -strand - -bg -ibam $bam > $negativeCoordinates0

	#1-BASED COORDINATES OF	rNMPs:
	#Obtain coordinates of rNMPs (3’ end of aligned read):
	bedtools genomecov -3 -strand + -d -ibam $bam > $positiveCoordinates1
	bedtools genomecov -3 -strand - -d -ibam $bam > $negativeCoordinates1

	#Remove rows where genome coverage equals 0
	awk '$3 != 0' $positiveCoordinates1 > temporary1 && mv temporary1 $positiveCoordinates1
	awk '$3 != 0' $negativeCoordinates1 > temporary2 && mv temporary2 $negativeCoordinates1

##########################################################################################################################################
	#STEP 3: Calculate background dNTP frequencies of reference genome

	#Location of input file
	referenceFasta1=$directory0/$subset.fa

	#Location of output directory
	output2=$directory/ribose-seq/results/Background-dNTP-Frequencies

	#Location of output file
	background=$output2/$reference.$subset.Background-dNTP-Frequencies.txt

	#Remove file if it already exists
	rm $background

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
	riboList=$output1/$sample.rNMP-list.$reference.$subset.txt
	riboFrequencies=$output1/$sample.rNMP-frequencies.$reference.$subset.txt

	#Remove files if they already exist
	rm $riboFrequencies $riboList

	#Select only rNMPs in subset:
	#Whole genome subset
	if [ $subset == "sacCer2" ] || [ $subset == "eColi" ] || [ $subset == "mm9" ] || [ $subset == "hg38" ] || [ $subset == "LL_1510A" ]; then
		awk -v "OFS=\t" '{print $4, $5}' $readCoordinates > temporary
	#Mitochondria subset
	elif [ $subset == "chrM" ]; then
    		grep 'chrM' $readCoordinates | awk -v "OFS=\t" '{print $4, $5}' - > temporary
	#Nuclear subset
	elif [ $subset == "nuclear" ]; then
		grep -v 'chrM' $readCoordinates | awk -v "OFS=\t" '{print $4, $5}' - > temporary
	fi

	#Print only rNMPs (3' end of reads):
	#rNMPs on positive strands (located at end of sequence)
	awk '$2 == "+" {print substr($0,length($0)-2)}' temporary > $riboList

	#rNMPs on negative strands (located at beginning of sequence)
	awk -v "OFS=\t" '$2 == "-" {print substr($0,0,1), $2}' temporary >> $riboList

	#Calculate count of each rNMP
	A_riboCount=$(awk '$1 == "A" && $2 == "+" || $1 == "T" && $2 == "-" {print $1, $2}' $riboList | wc -l)
	C_riboCount=$(awk '$1 == "C" && $2 == "+" || $1 == "G" && $2 == "-" {print $1, $2}' $riboList | wc -l)
	G_riboCount=$(awk '$1 == "G" && $2 == "+" || $1 == "C" && $2 == "-" {print $1, $2}' $riboList | wc -l)
	U_riboCount=$(awk '$1 == "T" && $2 == "+" || $1 == "A" && $2 == "-" {print $1, $2}' $riboList | wc -l)

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

	#Remove temporary file
	rm temporary

##########################################################################################################################################
	#STEP 5: Obtain coordinates and sequences of +/- 100 downstream/upstream dNTPs from rNMPs

	#Location of input files
	referenceBED=$directory0/$reference.bed
	referenceFasta2=$directory0/$reference.fa

	#Location of output directory
	output3=$directory2/dNTPs/$subset

	#Create directory if it does not already exist
	if [[ ! -d $output3 ]]; then
    		mkdir -p $output3
	fi	

	#Location of output files
	coordinates=$output1/$sample.rNMP-coordinates.bed
	
	positiveCoordinates=$output1/$sample.rNMP-coordinates.positive.bed
	negativeCoordinates=$output1/$sample.rNMP-coordinates.negative.bed
	
	positiveUpstreamIntervals=$output3/$sample.upstream-intervals.positive.bed
	negativeUpstreamIntervals=$output3/$sample.upstream-intervals.negative.bed
	positiveDownstreamIntervals=$output3/$sample.downstream-intervals.positive.bed
	negativeDownstreamIntervals=$output3/$sample.downstream-intervals.negative.bed
	
	positiveUpstreamSequences=$output3/$sample.upstream-sequences.positive.fa
	negativeUpstreamSequences=$output3/$sample.upstream-sequences.negative.fa
	positiveDownstreamSequences=$output3/$sample.downstream-sequences.positive.fa
	negativeDownstreamSequences=$output3/$sample.downstream-sequences.negative.fa
		
	#Obtain positions of rNMPs (3’ end of aligned reads)
	bedtools genomecov -3 -bg -ibam $bam > temporary1

	#Print only columns containing coordinates (eliminate column containing coverage values)
	awk -v "OFS=\t" '{print $1, $2, $3}' temporary1 > temporary2

	#Combine rNMP coordiantes with strand information in tabular format (i.e., 2micron 94 95 +)
	paste temporary2 $readCoordinates | awk -v "OFS=\t" '{print $1, $2, $3, $8}' > $coordinates
	
	awk -v "OFS=\t" '$4 == "+" {print $1, $2, $3, $4}' $coordinates > $positiveCoordinates
	awk -v "OFS=\t" '$4 == "-" {print $1, $2, $3, $4}' $coordinates > $negativeCoordinates
	
	#Obtain coordinates of dNTPs located 100 bp upstream from rNMPs:
	#Note: For positive strands, left = upstream. For negative strands, right = upstream
	bedtools flank -i $positiveCoordinates -g $referenceBED -l 2 -r 0 > $positiveUpstreamIntervals
	bedtools flank -i $negativeCoordinates -g $referenceBED -l 0 -r 2 > $negativeUpstreamIntervals
	
	#Obtain coordinates of dNTPs located 100 bp downstream from rNMPs:
	#Note: For positive strands, right = downstream. For negative strands, left = downstream
	bedtools flank -i $positiveCoordinates -g $referenceBED -l 0 -r 2 > $positiveDownstreamIntervals
	bedtools flank -i $negativeCoordinates -g $referenceBED -l 2 -r 0 > $negativeDownstreamIntervals

	#Obtain sequences of dNTPs located 100 bp upstream from rNMPs:
	bedtools getfasta -fi $referenceFasta2 -bed $positiveUpstreamIntervals -fo $positiveUpstreamSequences
	bedtools getfasta -fi $referenceFasta2 -bed $negativeUpstreamIntervals -fo $negativeUpstreamSequences
	
	#Obtain sequences of dNTPs located 100 bp downstream from rNMPs:
	bedtools getfasta -fi $referenceFasta2 -bed $positiveDownstreamIntervals -fo $positiveDownstreamSequences
	bedtools getfasta -fi $referenceFasta2 -bed $negativeDownstreamIntervals -fo $negativeDownstreamSequences

	rm temporary1 temporary2
##########################################################################################################################################
	#STEP 6: Tabulate sequences of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	#Location of output directory
	output4=$directory2/dNTPs/$subset/Columns/upstream
	output5=$directory2/dNTPs/$subset/Columns/downstream
	output6=$directory2/dNTPs/$subset/Columns/upstream/sequences
	output7=$directory2/dNTPs/$subset/Columns/downstream/sequences
				
	#Create directories if they do not already exist
	#if [[ ! -d $output4 && $output5 && $output6 && $output7 ]]; then
    	mkdir -p $output4 $output5 $output6 $output7
	#fi
	
	#Select all reads located in genomic DNA
	if [ $subset == "sacCer2" ] || [ $subset == "eColi" ] || [ $subset == "mm9" ] || [ $subset == "hg38" ] || [ $subset == "LL_1510A" ]; then
		cat $positiveUpstreamSequences > $output3/temporary1.positive.upstream
		cat $negativeUpstreamSequences > $output3/temporary1.negative.upstream
		cat $positiveDownstreamSequences > $output3/temporary1.positive.downstream
		cat $negativeDownstreamSequences > $output3/temporary1.negative.downstream
	#Select only reads located in nuclear DNA
	elif [ $subset == "nuclear" ]; then
		(grep -A 1 '>2micron' $positiveUpstreamSequences && grep -P -A 1 '>chr(?!M)' $positiveUpstreamSequences) > $output3/temporary1.positive.upstream
		(grep -A 1 '>2micron' $negativeUpstreamSequences && grep -P -A 1 '>chr(?!M)' $positiveUpstreamSequences) > $output3/temporary1.negative.upstream
		(grep -A 1 '>2micron' $positiveDownstreamSequences && grep -P -A 1 '>chr(?!M)' $positiveUpstreamSequences) > $output3/temporary1.positive.downstream
		(grep -A 1 '>2micron' $negativeDownstreamSequences && grep -P -A 1 '>chr(?!M)' $positiveUpstreamSequences) > $output3/temporary1.negative.downstream
	#Select only reads located in mitochondrial DNA
	elif [ $subset == "chrM" ]; then
		grep -A 1 '>chrM' $positiveUpstreamSequences > $output3/temporary1.positive.upstream
		grep -A 1 '>chrM' $negativeUpstreamSequences > $output3/temporary1.negative.upstream
		grep -A 1 '>chrM' $positiveDownstreamSequences > $output3/temporary1.positive.downstream
		grep -A 1 '>chrM' $negativeDownstreamSequences > $output3/temporary1.negative.downstream
	fi
				
	#Reverse complement upstream/downstream sequences on negative strand
	seqtk seq -r temporary1.negative.upstream > $output3/temporary2.negative.upstream
	seqtk seq -r temporary1.negative.downstream > $output3/temporary2.negative.downstream

	#Output only sequences in FASTA file (exclude all header lines)
	grep -v '>' temporary1.positive.upstream > $output3/temporary2.positive.upstream
	grep -v '>' temporary2.negative.upstream > $output3/temporary3.negative.upstream
	grep -v '>' temporary1.positive.downstream > $output3/temporary2.positive.downstream
	grep -v '>' temporary2.negative.downstream > $output3/temporary3.negative.downstream
				
	#Reverse upstream/downstream sequences on negative strand
	cat temporary3.negative.upstream|rev > $output3/temporary4.negative.upstream
	cat temporary3.negative.downstream|rev > $output3/temporary4.negative.downstream
				
	sequences1=$output6/$sample.upstream-sequences.$reference.$subset.txt
	sequences2=$output7/$sample.downstream-sequences.$reference.$subset.txt
	
	rm $sequences1 $sequences2
				
	#Combine upstream (+/-) and downstream (+/-) sequences into two files
	cat temporary2.positive.upstream temporary4.negative.upstream > $sequences1
	cat temporary2.positive.downstream temporary4.negative.downstream > $sequences2
			
	columns1=$output6/$sample.upstream-sequences.$reference.$subset.tab
	columns2=$output7/$sample.downstream-sequences.$reference.$subset.tab
				
	#Insert tabs between each nucleotide
	cat $sequences1 | sed 's/.../& /2g;s/./& /g' > $columns1
	cat $sequences2 | sed 's/.../& /2g;s/./& /g' > $columns2

	for i in {1..2}; do
		#Location of output files
		lists1=$output4/$sample.column.$i.upstream.$reference.$subset.txt
		lists2=$output5/$sample.column.$i.downstream.$reference.$subset.txt
		#Save lists of dNTPs at each +/- 100 bp downstream/upstream position
		awk -v field=$i '{ print $field }' $columns1 > $lists1
		awk -v field=$i '{ print $field }' $columns2 > $lists2
	done

##########################################################################################################################################

	#STEP 7: Calculate frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs

	locations="upstream downstream"
	
	#Location of output directory
	output5=$directory2/dNTPs/$subset/Raw-Data

	#Create directory if it does not already exist
	mkdir -p $output5
		
	#Remove old files if they already exist
	rm $output5/*.txt

	for location in ${locations[@]}; do

		#Location of output files (indivdiual base frequencies)
		A_baseFrequencies=$output5/A_dNTP-frequencies.$reference.$subset.$location.txt
		C_baseFrequencies=$output5/C_dNTP-frequencies.$reference.$subset.$location.txt
		G_baseFrequencies=$output5/G_dNTP-frequencies.$reference.$subset.$location.txt
		T_baseFrequencies=$output5/T_dNTP-frequencies.$reference.$subset.$location.txt

		#Location of output file (combined base frequencies)
		baseFrequencies=$output5/$sample.dNTP-frequencies.$reference.$subset.$location.txt
	
		#Calculate dNTP frequencies for each +/- 100 downstream/upstream position
		for file in $directory2/dNTPs/$subset/Columns/$location/$sample*.txt; do

			#Calculate count of each dNTP
			A_baseCount=$(grep -v '>' $file | grep -o 'A' - | wc -l)
			C_baseCount=$(grep -v '>' $file | grep -o 'C' - | wc -l)
			G_baseCount=$(grep -v '>' $file | grep -o 'G' - | wc -l)
			T_baseCount=$(grep -v '>' $file | grep -o 'T' - | wc -l)

			#Calculate total number of dNTPs
			total_baseCount=$(($A_baseCount+$C_baseCount+$G_baseCount+$T_baseCount))
	
			#Calculate raw frequencies of dNTPs
			A_rawBaseFrequency=$(echo "scale = 12; $A_baseCount/$total_baseCount" | bc | awk '{printf "%.12f\n", $0}')
			C_rawBaseFrequency=$(echo "scale = 12; $C_baseCount/$total_baseCount" | bc | awk '{printf "%.12f\n", $0}')
			G_rawBaseFrequency=$(echo "scale = 12; $G_baseCount/$total_baseCount" | bc | awk '{printf "%.12f\n", $0}')
			T_rawBaseFrequency=$(echo "scale = 12; $T_baseCount/$total_baseCount" | bc | awk '{printf "%.12f\n", $0}')

			#Calculate normalized frequencies of dNTPs
			A_baseFrequency=$(echo "scale = 12; $A_rawBaseFrequency/$A_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
			C_baseFrequency=$(echo "scale = 12; $C_rawBaseFrequency/$C_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
			G_baseFrequency=$(echo "scale = 12; $G_rawBaseFrequency/$G_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
			T_baseFrequency=$(echo "scale = 12; $T_rawBaseFrequency/$T_backgroundFrequency" | bc | awk '{printf "%.12f\n", $0}')
		
			#Save normalized frequencies of dNTPs to TXT file
			echo $A_baseFrequency >> $A_baseFrequencies
			echo $C_baseFrequency >> $C_baseFrequencies
			echo $G_baseFrequency >> $G_baseFrequencies
			echo $T_baseFrequency >> $T_baseFrequencies

			#Remove old file if it already exists
			if [ -e "$baseFrequencies" ]; then
    				rm $baseFrequencies
			fi

			#Save frequencies of dNTPs located +/- 100 base pairs downstream/upstream from rNMPs to one TXT file
			paste $A_baseFrequencies $C_baseFrequencies $G_baseFrequencies $T_baseFrequencies >> $baseFrequencies
		done
	done

##########################################################################################################################################
	#STEP 8: Create dataset file containing nucleotide frequencies needed for plotting

	#Location of input files
	upstreamBaseFrequencies=$output5/$sample.dNTP-frequencies.$reference.$subset.upstream.txt
	downstreamBaseFrequencies=$output5/$sample.dNTP-frequencies.$reference.$subset.downstream.txt
	
	#Location of output directory
	output6=$directory2/Datasets/$subset

	#Create directory if it does not already exist
    	mkdir -p $output6

	#Location of output file
	dataset=$output6/$sample.nucleotide-frequencies-dataset.$reference.$subset.txt
	zoomed=$output6/$sample.nucleotide-frequencies-zoomed.$reference.$subset.txt
	
	#Remove old file if it already exists
	rm $dataset

	#Print values -100 to 100
	#seq -100 1 100 > temporary1
	seq -2 1 2 > temporary1
	
	#Save files containing rNMP and upstream/downstream dNTP frequencies to one file
	cat $upstreamBaseFrequencies $riboFrequencies $downstreamBaseFrequencies >> temporary2

	#Save files to one combined TXT file
	paste temporary1 temporary2 > temporary3

	#Add header line containing nucleotides to beginning of file 
	echo -e "\tA\tC\tG\tU/T" > $dataset; cat temporary3 >> $dataset;
	
	#Smaller dataset
	head -117 $dataset | tail -31 > temporary4
	echo -e "\tA\t\tC\t\tG\t\tU/T" > $zoomed; cat temporary4 >> $zoomed;

	#Remove temporary files
	rm temporary1 temporary2 temporary3 temporary4

	echo "Calculation of nucleotide frequencies for $sample is complete"

done
