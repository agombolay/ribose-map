#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates the ribonucleotide frequencies located at 3' position of input BED file

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

#Location of "Reference" directory
directory0=$directory/ribose-seq/reference/

#Location of "Alignment" directory
directory1=$directory/ribose-seq/results/$reference/$sample/Alignment

#Location of "dNTP-Frequencies" directory
directory2=$directory/ribose-seq/results/$reference/$sample/dNTP-Frequencies

##############################################################################################################################
#STEP 1: Covert BAM alignment file to FASTA format

#Location of input file
bam=$directory1/$sample.bam

#Location of output directory
output1=$directory2/Ribonucleotides/$subset

#Create directory for output if it does not already exist
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

##############################################################################################################################
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
#Obtain positions of rNMPs (3’ end of aligned read) for positive strand:
bedtools genomecov -3 -strand + -bg -ibam $bam > $positiveCoordinates0

#Obtain positions of rNMPs (3’ end of aligned read) for negative strand:
bedtools genomecov -3 -strand - -bg -ibam $bam > $negativeCoordinates0

#1-BASED COORDINATES OF	rNMPs:
#Obtain positions of rNMPs (3’ end of aligned read) for positive strand:
bedtools genomecov -3 -strand + -d -ibam $bam > $positiveCoordinates1

#Obtain positions of rNMPs (3’ end of aligned read) for negative strand:
bedtools genomecov -3 -strand - -d -ibam $bam > $negativeCoordinates1

#Remove rows where genome coverage equals 0
awk '$3 != 0' $positiveCoordinates1 > temporary1.txt

#Remove rows where genome coverage equals 0
awk '$3 != 0' $negativeCoordinates1 > temporary2.txt

#Change filename back to original
mv temporary1.txt $positiveCoordinates1

#Change filename back to original
mv temporary2.txt $negativeCoordinates1

##############################################################################################################################
#STEP 3: Calculate background dNTP frequencies of reference genome

#Location of input file
referenceFasta=$directory0/$subset.fa

#Location of output directory
output2=$directory/ribose-seq/results/Background-dNTP-Frequencies

#Location of output file
background=$output2/$reference.$subset.Background-dNTP-Frequencies.txt

#Remove file if it already exists
rm $background

#Calculate counts of each nucleotide
A_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'A' - | wc -l)
C_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'C' - | wc -l)
G_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'G' - | wc -l)
T_backgroundCount=$(grep -v '>' $referenceFasta | grep -o 'T' - | wc -l)

#Calculate total number of nucleotides
total_backgroundCount=$(($A_backgroundCount+$C_backgroundCount+$G_backgroundCount+$T_backgroundCount))

#Calculate frequency of each nucleotide
A_backgroundFrequency=$(bc <<< "scale = 4; `expr $A_backgroundCount/$total_backgroundCount`")
C_backgroundFrequency=$(bc <<< "scale = 4; `expr $C_backgroundCount/$total_backgroundCount`")
G_backgroundFrequency=$(bc <<< "scale = 4; `expr $G_backgroundCount/$total_backgroundCount`")
T_backgroundFrequency=$(bc <<< "scale = 4; `expr $T_backgroundCount/$total_backgroundCount`")

#Print nucleotide frequencies to output file
echo "A Background Frequency: $A_backgroundFrequency" >> $background
echo "C Background Frequency: $C_backgroundFrequency" >> $background
echo "G Background Frequency: $G_backgroundFrequency" >> $background
echo "T Background Frequency: $T_backgroundFrequency" >> $background

##############################################################################################################################
#STEP 4: Calculate Ribonucleotide Frequencies

#Location of output files
list=$output1/$sample.rNMP-list.$reference.$subset.txt
frequencies=$output1/$sample.rNMP-frequencies.$reference.$subset.txt

#Remove file if it already exists
rm $frequencies $list

#Select only ribonucleotides of genome subset:
#Whole genome subset
if [[ $subset == "sacCer2" ]]; then
	awk -v "OFS=\t" '{print $4, $5}' $readCoordinates > temporary.txt
#Mitochondria subset
elif [[ $subset == "mitochondria" ]]; then
    	grep 'chrM' $readCoordinates | awk -v "OFS=\t" '{print $4, $5}' - > temporary.txt
#Nuclear subset
elif [[ $subset == "nuclear" ]]; then
	grep -v 'chrM' $readCoordinates | awk -v "OFS=\t" '{print $4, $5}' - > temporary.txt
fi

#Print only ribonucleotides (3' end of read):
#Ribonucleotides on positive strands (located at end of sequence)
awk '$2 == "+" {print substr($0,length($0)-2)}' temporary.txt > $list

#Ribonucleotides on negative strands (located at beginning of sequence)
awk -v "OFS=\t" '$2 == "-" {print substr($0,0,1), $2}' temporary.txt >> $list

#Calculate count of each ribonucleotide
A_riboCount=$(awk '$1 == "A" && $2 == "+" || $1 == "T" && $2 == "-" {print $1, $2}' $list | wc -l)
C_riboCount=$(awk '$1 == "C" && $2 == "+" || $1 == "G" && $2 == "-" {print $1, $2}' $list | wc -l)
G_riboCount=$(awk '$1 == "G" && $2 == "+" || $1 == "C" && $2 == "-" {print $1, $2}' $list | wc -l)
U_riboCount=$(awk '$1 == "T" && $2 == "+" || $1 == "A" && $2 == "-" {print $1, $2}' $list | wc -l)

#Calculate total number of ribonucleotides
total_riboCount=$(($A_riboCount+$C_riboCount+$G_riboCount+$U_riboCount))

#Calculate raw frequency of each ribonucleotide
A_rawFrequency=$(bc <<< "scale = 4; `expr $A_riboCount/$total_riboCount`")
C_rawFrequency=$(bc <<< "scale = 4; `expr $C_riboCount/$total_riboCount`")
G_rawFrequency=$(bc <<< "scale = 4; `expr $G_riboCount/$total_riboCount`")
U_rawFrequency=$(bc <<< "scale = 4; `expr $U_riboCount/$total_riboCount`")

#Calculate normalized frequency of each ribonucleotide
A_riboFrequency=$(bc <<< "scale = 4; `expr $A_rawFrequency/$A_backgroundFrequency`")
C_riboFrequency=$(bc <<< "scale = 4; `expr $C_rawFrequency/$C_backgroundFrequency`")
G_riboFrequency=$(bc <<< "scale = 4; `expr $G_rawFrequency/$G_backgroundFrequency`")
U_riboFrequency=$(bc <<< "scale = 4; `expr $U_rawFrequency/$T_backgroundFrequency`")

#Print ribonucleotide frequencies to output file
echo "$A_riboFrequency $C_riboFrequency $G_riboFrequency $U_riboFrequency" | column  > $frequencies

#Remove temporary file
rm temporary.txt

##############################################################################################################################
#STEP 5: Obtain coordinates and sequences of +/- 100 downstream/upstream dNTPs from rNMPs

#Location of input file
referenceBED=$directory0/$reference.bed

#Location of output directory
output3=$directory2/Nucleotides/$subset

#Create directory for output if it does not already exist
if [[ ! -d $output3 ]]; then
    	mkdir -p $output3
fi	

#Location of output files
coordinates=$output1/$sample.rNMP-coordinates.bed
upstreamIntervals=$output3/$sample.upstream-intervals.bed
upstreamSequences=$output3/$sample.upstream-sequences.tab
downstreamIntervals=$output3/$sample.downstream-intervals.bed
downstreamSequences=$output3/$sample.downstream-sequences.tab

temporary1=$output3/temporary1.bed
temporary2=$output3/temporary2.bed

#Obtain positions of rNMPs (3’ end of aligned reads)
bedtools genomecov -3 -bg -ibam $bam > $coordinates

#Remove column containing coverage values
awk '!($4="")' $coordinates > $temporary1

#Then, change file back to its original name
mv $temporary1 $coordinates

#Make sure file is tab-delimited
sed 's/ \+/\t/g' $coordinates > $temporary2

#Then, change file back to its original name
mv $temporary2 $coordinates

#Obtain coordinates of +/- 100 downstream/upstream dNTPs from rNMPs:
bedtools flank -i $coordinates -g $referenceBED -l 100 -r 0 > $upstreamIntervals
bedtools flank -i $coordinates -g $referenceBED -l 0 -r 100 > $downstreamIntervals

#Obtain sequences of +/- 100 downstream/upstream dNTPs from rNMPs:
bedtools getfasta -fi $referenceFasta -bed $upstreamIntervals -tab -fo $upstreamSequences
bedtools getfasta -fi $referenceFasta -bed $downstreamIntervals -tab -fo $downstreamSequences

##############################################################################################################################
#STEP 6: Output sequences of +/- 100 downstream/upstream dNTPs from rNMPs into tabular format

locations="upstream downstream"

for location in ${locations[@]}; do

	for file in "$output3/$sample.$location.sequences.tab"; do
		
		#Location of output directory
		output4=$directory2/Nucleotides/$subset/Columns/$location

		#Create directory for output if it does not already exist
		if [[ ! -d $output4 ]]; then
    			mkdir -p $output4
		fi

		#Create directory for output if it does not already exist
		if [[ ! -d $output4/sequences ]]; then
                	mkdir -p $output4/sequences
        	fi

		#Location of output files
		selection=$output4/sequences/$sample.$location-sequences.$subset.txt
		sequences=$output4/sequences/$sample.$location-sequences.$subset.raw.txt
		columns=$output4/sequences/$sample.$location-sequences.$subset.columns.txt

		if [ $subset == "sacCer2" ]; then
			cat $file > $selection
		elif [ $subset == "mitochondria" ]; then
			grep 'chrM' $file > $selection
		elif [ $subset == "nuclear" ]; then
			grep -v 'chrM' $file > $selection
		fi

		#Print sequences to new file
		awk -v "OFS=\t" '{print $2}' $selection > $sequences

		#Insert tabs between each nucleotide
		cat $sequences | sed 's/.../& /2g;s/./& /g' > $columns

		for i in {1..100}; do
			awk -v field=$i '{ print $field }' $columns > $output4/$sample.column.$i.$location.$subset.txt
		done
	done
done

##############################################################################################################################
#STEP 7: Calculate frequencies of +/- 100 downstream/upstream dNTPs from rNMPs

#Location of output directory
output5=$directory2/Nucleotides/$subset/Raw-Data

#Create directory for output if it does not already exist
if [[ ! -d $output5 ]]; then
    	mkdir -p $output5
fi
		
#Remove old .txt files
rm $output5/*.txt

for location in ${locations[@]}; do

	A_dNTP_frequencies=$output5/A_dNTP-frequencies.$reference.$subset.$location.txt
	C_dNTP_frequencies=$output5/C_dNTP-frequencies.$reference.$subset.$location.txt
	G_dNTP_frequencies=$output5/G_dNTP-frequencies.$reference.$subset.$location.txt
	T_dNTP_frequencies=$output5/T_dNTP-frequencies.$reference.$subset.$location.txt
	
	dNTP_Frequencies=$output5/$sample.dNTP-frequencies.$reference.$subset.$location.txt
		
	input=$directory2/Nucleotides/$subset/Columns/$location/$sample*.txt
	
	for file in ${input[@]}; do
	
		A_dNTP_count=$(grep -v '>' $file | grep -o 'A' - | wc -l)
		C_dNTP_count=$(grep -v '>' $file | grep -o 'C' - | wc -l)
		G_dNTP_count=$(grep -v '>' $file | grep -o 'G' - | wc -l)
		T_dNTP_count=$(grep -v '>' $file | grep -o 'T' - | wc -l)

		total_dNTP_count=$(($A_dNTP_count+$C_dNTP_count+$G_dNTP_count+$T_dNTP_count))
	
		A_dNTP_rawFrequency=$(bc <<< "scale = 4; `expr $A_dNTP_count/$total_dNTP_count`")
		C_dNTP_rawFrequency=$(bc <<< "scale = 4; `expr $C_dNTP_count/$total_dNTP_count`")
		G_dNTP_rawFrequency=$(bc <<< "scale = 4; `expr $G_dNTP_count/$total_dNTP_count`")
		T_dNTP_rawFrequency=$(bc <<< "scale = 4; `expr $T_dNTP_count/$total_dNTP_count`")

		A_dNTP_frequency=$(bc <<< "scale = 4; `expr $A_dNTP_rawFrequency/$A_backgroundFrequency`")
        	C_dNTP_frequency=$(bc <<< "scale = 4; `expr $C_dNTP_rawFrequency/$C_backgroundFrequency`")
        	G_dNTP_frequency=$(bc <<< "scale = 4; `expr $G_dNTP_rawFrequency/$G_backgroundFrequency`")
        	T_dNTP_frequency=$(bc <<< "scale = 4; `expr $T_dNTP_rawFrequency/$T_backgroundFrequency`")

		echo $A_dNTP_frequency >> $A_dNTP_frequencies
		echo $C_dNTP_frequency >> $C_dNTP_frequencies
		echo $G_dNTP_frequency >> $G_dNTP_frequencies
		echo $T_dNTP_frequency >> $T_dNTP_frequencies

		if [ -e "$dNTP_Frequencies" ]; then
    			rm $dNTP_Frequencies
		fi

		paste $A_dNTP_frequencies $C_dNTP_frequencies $G_dNTP_frequencies $T_dNTP_frequencies >> $dNTP_Frequencies
	done
done

##############################################################################################################################
#STEP 8: Create dataset file containing nucleotide frequencies needed for plotting

#Location of output directory
output6=$directory2/Datasets/$subset

#Location of output file
dataset=$output6/$sample.nucleotide-frequencies-dataset.$subset.txt

#Create directory for output if it does not already exist
if [[ ! -d $output6 ]]; then
    	mkdir -p $output6
fi

#Remove old .txt files
rm $output6/*.txt

#Print values -100 to 100
seq -100 1 100 > temporary1.txt

#Combine nucleotide frequency files
cat $output5/$sample.Normalized_Nucleotide_Frequencies.$subset.upstream.txt \
$output1/$sample.$reference.$subset.ribonucleotide-frequencies.txt \
$output5/$sample.Normalized_Nucleotide_Frequencies.$subset.downstream.txt >> temporary2.txt

#Merge two files into final .txt file
paste temporary1.txt temporary2.txt > temporary3.txt

#Add Header to beginning of .txt file 
echo "Position A C G U/T" | awk '{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5}' \
| cat - temporary3.txt > temp && mv temp temporary3.txt

#Make sure file is tab-delimited
column -t temporary3.txt > $dataset

#Remove temporary files
rm temporary1.txt temporary2.txt temporary3.txt
