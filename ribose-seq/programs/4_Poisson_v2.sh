#!/usr/bin/env bash

#Author: Alli Gombolay
#This program counts the number of positions in the genome with 0...X number of rNMPs

#COMMAND LINE OPTIONS

#Usage statement of the program
function usage () {
	echo "Usage: 4_Poisson.sh [-i] 'Sample' [-r] 'Reference' [-s] 'Subset' [-d] 'Directory' [-h]
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

#############################################################################################################################

#Input files
bed=$directory/ribose-seq/reference/$reference.bed
bam=$directory/ribose-seq/results/$reference/$sample/Alignment/$sample.bam

#Output files
coverage=$directory/ribose-seq/results/$reference/$sample/Poisson/$sample.rNMP-coverage.bed
counts1=$directory/ribose-seq/results/$reference/$sample/Poisson/$sample.rNMP-counts.bed

#Obtain coverage at 3' positions of BAM file
if [ $subset == "nuclear" ]; then
	#Calculate total number of genome positions
	total=$(grep -v 'chrM' $bed | awk '{sum+=$2} END{print sum}' -)
	#Select only nuclear DNA regions from BED file
	bedtools genomecov -3 -dz -ibam $bam -g $bed | grep -v 'chrM' - > $coverage
elif [ $subset == "chrM" ]; then
	#Calculate total number of genome positions
	total=$(grep 'chrM' $bed | awk '{print $2}' -)
	#Select only mitochondrial DNA regions from BED file
	bedtools genomecov -3 -dz -ibam $bam -g $bed | grep 'chrM' - > $coverage
fi

#Maximum value of genome coverage in BED file
#maximum=$(sort -nk 3 $coverage | tail -1 - | awk '{print $3}' -)

#Number of positions with X number of rNMPs
#for i in $(seq 1 $maximum); do
#	positions1+=($(awk '$3 == ('$i')' $coverage | wc -l))
#done

#Number of positions with 0 rNMPs (total positions-positions with rNMPs)
#positions2=$(echo "($total-$(wc -l $coverage | awk '{print $1}' -))" | bc)

#Print observed count data to fit to Poisson distribution
#( IFS=$'\n'; echo -e "$positions2\n${positions1[*]}" ) > $counts1

#echo -e "rNMPs\tPositions"
#echo -e "0\t$positions2"

#array1=($(seq 1 $maximum))
#array2=(${positions1[*]})

#sum=0
#for i in "${!array1[@]}"; do
#        echo -e "${array1[i]}\t${array2[i]}"
#        (( sum +=$(echo "${array1[$i]}*${array2[$i]}" | bc) ))
#done

#Calculate lambda for Poisson Distribution (in scientific notation)
#lambda=$(echo "scale = 12; $sum/$total" | bc | awk '{printf "%e\n", $0}')

#echo "Total rNMPs:" $sum
#echo "Lambda:" $lambda

#Version 2: Proportion of windows that have x number of ribos

#Input files
bed=$directory/ribose-seq/reference/$reference.bed
sorted=$directory/ribose-seq/results/$reference/$sample/Coordinates/$subset/$sample.rNMP-coordinates.sorted.bed

#Output directories
output1=$directory/ribose-seq/reference/
output2=$directory/ribose-seq/results/$reference/$sample/Poisson

#Create directories
mkdir -p $output1 $output2

#Output files
binnedData=$output2/$sample.binned.data.bed
referenceWindows=$output1/$reference.windows.bed

counts1=$output2/$sample.probability-mass-function.counts.txt
counts2=$output2/$sample.cumulative-distribution.counts.txt

proportions1=$output2/$sample.probability-mass-function.proportions.txt
proportions2=$output2/$sample.cumulative-distribution.proportions.txt

#Separate reference genome into 2.5 kb windows
bedtools makewindows -g $bed -w 2500 > $referenceWindows

#Select only data of interest
if [ $subset == "nuclear" ]; then
	#Determine regions of BED files that intersect and count number of overlaps (nuclear)
	bedtools intersect -a $referenceWindows -b $sorted -c -sorted -nonamecheck | grep -v 'chrM' - > $binnedData
elif [ $subset == "chrM" ]; then
	#Determine regions of BED files that intersect and count number of overlaps (chrM)
	bedtools intersect -a $referenceWindows -b $sorted -c -sorted -nonamecheck | grep 'chrM' - > $binnedData
fi

#Maximum value of genome coverage in BED file
maximum=$(sort -nk 4 $binnedData | tail -1 - | awk '{print $4}' -)

#Number of positions with X number of rNMPs
for i in $(seq 0 $maximum); do
	windows+=($(awk '$4 == ('$i')' $binnedData | wc -l))
done

#Probability Mass Function Counts
#Print number of windows in genomic region with exactly 0...X rNMPs
paste <(echo "$(seq 0 $maximum)") <(cat <( IFS=$'\n'; echo "${windows[*]}" )) > $counts1

#Cumulative Distribution Counts
#Print number of windows in genomic region with greater than or equal to 0...X rNMPs
for i in $(seq $(wc -l < $counts1) -1 1); do
	head -$(wc -l < $counts1) $counts1 | tail -${i} | awk '{ SUM += $2} END { print SUM }' >> $counts2
done

#Total number of windows
total=$(awk '{ SUM += $2} END { print SUM }' $counts1)

#Proportions of windows (P(X=x))
for i in ${windows[*]}; do
	values1+=($(echo "scale = 12; ($i/$total)" | bc | awk '{printf "%.12f\n", $0}'))
done

#Proportions of windows (P(X>=x))
for i in ${values1[*]}; do
	values2+=($(echo "scale = 12; (1-$i)" | bc | awk '{printf "%.12f\n", $0}'))
done

#Print proportions of windows (probability mass function and cumulative distribution)
paste <(echo "$(seq 0 $maximum)") <(cat <( IFS=$'\n'; echo "${values1[*]}" )) > $proportions1
paste <(echo "$(seq 0 $maximum)") <(cat <( IFS=$'\n'; echo "${values2[*]}" )) > $proportions2

#Calculate lambda for Poisson Distribution
echo "Lambda:" $(echo "scale = 12; $(wc -l < $sorted)/$total" | bc | awk '{printf "%.5f\n", $0}')
#echo "Lambda:" $lambda

#variable=0
#proportions=()
#counts0=$(awk '$4 == 0' FS15.trimmed.v1.binned.data.bed | wc -l)
#counts2=$(awk '$4 > 9' FS15.trimmed.v1.binned.data.bed | awk '{sum+=$4} END{print sum}')

#for i in {1..9}; do
#	(( variable+=$(awk '$4 == ('$i')' FS15.trimmed.v1.binned.data.bed | awk '{sum+=$4} END{print sum}') ))
#	counts1+=($(awk '$4 == ('$i')' FS15.trimmed.v1.binned.data.bed | awk '{sum+=$4} END{print sum}'))
#done

#total=$(($counts0+$variable+$counts2))
#for value in $counts0 ${counts1[*]}; do
#	proportions+=($(echo "scale = 12; ($value/$total)" | bc | awk '{printf "%.12f\n", $0}'))
#done

#for value in ${proportions[*]}; do
#	final+=($(echo "scale = 12; (1-$value)" | bc | awk '{printf "%.12f\n", $0}'))
#done
