#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Creates bedgraph files for forward and reverse strands
#2. Saves coverage of rNMPs per chromosome to separate files

. /data2/users/agombolay3/Ribose-Map/config.txt

output=$directory/results/$sample/distribution
mkdir -p $output

#############################################################################################################################
#Create BED file for reference genome
cut -f 1,2 $reference.fai > $output/reference.bed
	
#Create file of rNMP coverage at chromosome coordinates
uniq -c $directory/results/$sample/coordinates/$sample.bed \
| awk -v "OFS=\t" '{print $2,$3,$4,$1,$5}' > $output/temp1.txt
	
#Save coverage of rNMPs per chromosome to separate files
for chromosome in $( awk '{print $1}' $output/reference.bed ); do
	grep -w "$chromosome" $output/temp1.txt > $output/$sample.$chromosome.bed
done
		
#Add trackline for forward strand to input into UCSC genome browser
echo "track type=bedGraph name="$sample-ForwardStrand" description="$sample-ForwardStrand" \
color=0,128,0 visibility=full" > $output/$sample-Forward.bg
		
#Add trackline for reverse strand to input into UCSC genome browser
echo "track type=bedGraph name="$sample-ReverseStrand" description="$sample-ReverseStrand" \
color=0,0,255 visibility=full" > $output/$sample-Reverse.bg
		
#Rearrange forward strand file so format is the same as bedgraph format
awk -v "OFS=\t" '$5 == "+" {print $0}' $output/temp1.txt >> $output/$sample-Forward.bg
#awk -v "OFS=\t" '$5 == "+" {print $2,$3,$4,$1}' $output/temp1.txt >> $output/$sample-Forward.bg

#Rearrange reverse strand file so format is the same as bedgraph format
awk -v "OFS=\t" '$5 == "-" {print $0}' $output/temp1.txt >> $output/$sample-Reverse.bg
#awk -v "OFS=\t" '$5 == "+" {print $2,$3,$4,$1}' $output/temp1.txt >> $output/$sample-Forward.bg

#############################################################################################################################
#Print status
echo "Status: Coverage of rNMPs have been determined for $sample"
	
#Remove temporary files
rm $output/reference.bed $output/temp1.txt
