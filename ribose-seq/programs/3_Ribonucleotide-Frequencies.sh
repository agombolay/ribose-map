#!/usr/bin/env bash
#Author: Alli Gombolay
#This program calculates the ribonucleotide frequencies located at 3' position of input BED file

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 5_Ribonucleotide-Frequencies.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
	-i Sample name (FS1, etc.)
	-s Subset of genome (sacCer2, nuclear, chrM)
	-r Reference genome assembly (sacCer2, etc.)
	-d Local directory ('/projects/home/agombolay3/data/repository/Ribose-seq-Project')"
}

#Use getopts function to create the command-line options ([-i], [-s], [-r], [-d], and [-h])
while getopts "i:s:r:d:h" opt;
do
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
if [ "$1" == "-h" ];
then
        exit
fi

##############################################################################################################################
#STEP 1: Covert BAM file to FASTA format

#Location of input BAM file
bam=$directory/ribose-seq/results/$reference/$sample/Alignment/$sample.bam

#Location of output directory
output=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Ribonucleotides

#Create directory for output if it does not already exist
if [[ ! -d $output ]];
then
	mkdir -p $output
fi
	
fastq=$output/$sample.fastq
fasta=$output/$sample.fasta

samtools bam2fq $bam > $fastq
seqtk seq -A $fastq > $fasta

##############################################################################################################################
#STEP 2: Obtain Ribonucleotide Coordinates

#Location of output files
bed=$output/$sample.bed
sam=$output/$sample.sam
coordinates=$output/$sample.coordinates.bed
positionsPositive0=$output/$sample.rNMPs.positive.0-based.txt
positionsNegative0=$output/$sample.rNMPs.negative.0-based.txt
positionsPositive1=$output/$sample.rNMPs.positive.1-based.txt
positionsNegative1=$output/$sample.rNMPs.negative.1-based.txt
	
#COORDINATES (0-BASED) of SEQUENCING READS

#Covert file from BAM to BED format using BEDtools software
bedtools bamtobed -i $bam > $bed

#Convert file from BAM to SAM format using SAMtools software
samtools view $bam > $sam

#Extract read coordinates, sequences, and strands from BED and SAM files and save it to new file
paste $bed $sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' > $coordinates

#0-BASED COORDINATES OF rNMPs

#Obtain positions of rNMPs (3’ end of each mapped read) for positive strand:
bedtools genomecov -3 -strand + -bg -ibam $bam > $positionsPositive0

#Obtain positions of rNMPs (3’ end of each mapped read) for negative strand:
bedtools genomecov -3 -strand - -bg -ibam $bam > $positionsNegative0

#1-BASED COORDINATES OF	rNMPs

#Obtain positions of rNMPs (3’ end of each mapped read) for positive strand:
bedtools genomecov -3 -strand + -d -ibam $bam > $positionsPositive1

#Remove rows where genome coverage equals 0
awk '$3 != 0' $positionsPositive1 > temporary

#Change filename back to original
mv temporary $positionsPositive1

#Obtain positions of rNMPs (3’ end of each mapped read) for negative strand:
bedtools genomecov -3 -strand - -d -ibam $bam > $positionsNegative1

#Remove rows where genome coverage equals 0
awk '$3 != 0' $positionsNegative1 > temporary

#Change filename back to original
mv temporary $positionsNegative1

##############################################################################################################################
#STEP 3: Calculate Ribonucleotide Frequencies

#Location of input BED file
bed=$directory/ribose-seq/results/$reference/$sample/Nucleotide-Frequencies/Ribonucleotides/$sample.coordinates.bed

#Print only ribonucleotides of genome subset to output file
#Whole genome subset
if [[ $subset == "sacCer2" ]];
then
	awk -v "OFS=\t" '{print $4, $5}' $bed > temporary.txt
#Mitochondria subset
elif [[ $subset == "mitochondria" ]];
then
    	grep 'chrM' $bed | awk -v "OFS=\t" '{print $4, $5}' - > temporary.txt
#Nuclear subset
elif [[ $subset == "nuclear" ]];
then
	grep -v 'chrM' $bed | awk -v "OFS=\t" '{print $4, $5}' - > temporary.txt
fi

#Print only ribonucleotides (3' end of read (end for + strand and start for - strand)) to output file

#Print ribonucleotides for positive strands (located at end of sequence)
awk '$2 == "+" {print substr($0,length($0)-2)}' temporary.txt > Ribonucleotide_List.$subset.txt

#Print ribonucleotides for negative strands (located at beginning of sequence)
awk -v "OFS=\t" '$2 == "-" {print substr($0,0,1), $2}' temporary.txt >> Ribonucleotide_List.$subset.txt

#Calculate count of "A" ribonucleotides
A_ribonucleotide_count=$(awk '$1 == "A" && $2 == "+" || $1 == "T" && $2 == "-" {print $1, $2}' Ribonucleotide_List.$subset.txt | wc -l)

#Calculate count of "C"	ribonucleotides
C_ribonucleotide_count=$(awk '$1 == "C" && $2 == "+" || $1 == "G" && $2 == "-" {print $1, $2}' Ribonucleotide_List.$subset.txt | wc -l)

#Calculate count of "G"	ribonucleotides
G_ribonucleotide_count=$(awk '$1 == "G" && $2 == "+" || $1 == "C" && $2 == "-" {print $1, $2}' Ribonucleotide_List.$subset.txt | wc -l)

#Calculate count of "U"	ribonucleotides
U_ribonucleotide_count=$(awk '$1 == "T" && $2 == "+" || $1 == "A" && $2 == "-" {print $1, $2}' Ribonucleotide_List.$subset.txt | wc -l)

total_ribonucleotide_count=$(($A_ribonucleotide_count+$C_ribonucleotide_count+$G_ribonucleotide_count+$U_ribonucleotide_count))

A_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $A_ribonucleotide_count/$total_ribonucleotide_count`")
C_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $C_ribonucleotide_count/$total_ribonucleotide_count`")
G_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $G_ribonucleotide_count/$total_ribonucleotide_count`")
U_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $U_ribonucleotide_count/$total_ribonucleotide_count`")
		
#A_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $A_ribonucleotide_frequency/$A_background_frequency`")
#C_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $C_ribonucleotide_frequency/$C_background_frequency`")
#G_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $G_ribonucleotide_frequency/$G_background_frequency`")
#U_normalized_ribonucleotide_frequency=$(bc <<< "scale = 4; `expr $U_ribonucleotide_frequency/$T_background_frequency`")

echo $A_ribonucleotide_count
echo $C_ribonucleotide_count
echo $G_ribonucleotide_count
echo $U_ribonucleotide_count

#echo $A_ribonucleotide_normalized_frequency
#echo $C_ribonucleotide_normalized_frequency
#echo $G_ribonucleotide_normalized_frequency
#echo $U_ribonucleotide_normalized_frequency

echo $total_ribonucleotide_count

rm temporary.txt
