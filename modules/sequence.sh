#!/usr/bin/env bash

#Author: Alli L. Gombolay
#Calculates frequencies of rNMP nucleotides and flanking nucleotides

#############################################################################################################################

#Load config file
. "$1"

#Create output directory and remove any old files
output=$repository/results/$sample/sequence$quality
rm -r $output; mkdir -p $output

#############################################################################################################################

#Calculate frequencies of rNMPs
for region in $unit "chromosomes"; do
		
	#Extract rNMP nucleotides from FASTA
	bedtools getfasta -s -fi $repository/$(basename $fasta .fa)-$region.fa -bed $repository/results/$sample/coordinate$quality/${sample}-$region.Combined.tab | grep -v '>' > $output/${sample}-$region.ribos.txt
			
	#Calculate counts of rNMPs
	A_Ribo=$(awk '$1 == "A" || $1 == "a"' $output/${sample}-$region.ribos.txt | wc -l)
	C_Ribo=$(awk '$1 == "C" || $1 == "c"' $output/${sample}-$region.ribos.txt | wc -l)
	G_Ribo=$(awk '$1 == "G" || $1 == "g"' $output/${sample}-$region.ribos.txt | wc -l)
	U_Ribo=$(awk '$1 == "T" || $1 == "t"' $output/${sample}-$region.ribos.txt | wc -l)
	
	#Calculate total number of rNMPs
	RiboTotal=$(($A_Ribo + $C_Ribo + $G_Ribo + $U_Ribo))
				
	#Calculate raw frequency of each rNMP
	A_RiboFreq=$(echo "($A_Ribo/$RiboTotal)" | bc -l)
	C_RiboFreq=$(echo "($C_Ribo/$RiboTotal)" | bc -l)
	G_RiboFreq=$(echo "($G_Ribo/$RiboTotal)" | bc -l)
	U_RiboFreq=$(echo "($U_Ribo/$RiboTotal)" | bc -l)
			
	for nuc in "A" "C" "G" "U" "Combined"; do

		#Save frequencies of rNMPs to TXT files
		if [[ $nuc == "A" ]]; then
			echo $A_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-Ribo$nuc.$region.$nuc.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.C.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.G.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.U.txt
				
			#Create BED file for only A ribonucleotide
			bedtools getfasta -s -fi $repository/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.Combined.tab | awk '$2 == "A" || $2 == "a"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-Ribo.$region.$nuc.bed

		elif [[ $nuc == "C" ]]; then
			echo $C_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-Ribo$nuc.$region.$nuc.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.A.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.G.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.U.txt
				
			#Create BED file for only C ribonucleotide
			bedtools getfasta -s -fi $repository/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.Combined.tab | awk '$2 == "C" || $2 == "c"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-Ribo.$region.$nuc.bed
				
		elif [[ $nuc == "G" ]]; then
			echo $G_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-Ribo$nuc.$region.$nuc.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.A.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.C.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.U.txt
				
			#Create BED file for only G ribonucleotide
			bedtools getfasta -s -fi $repository/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.Combined.tab | awk '$2 == "G" || $2 == "g"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-Ribo.$region.$nuc.bed

		elif [[ $nuc == "U" ]]; then
			echo $U_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-Ribo$nuc.$region.$nuc.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.A.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.C.txt
			echo 'NA' > $output/${sample}-Ribo$nuc.$region.G.txt
				
			#Create BED file for only U ribonucleotide
			bedtools getfasta -s -fi $repository/$(basename $fasta .fa)-$region.fa -tab -bed $repository/results/$sample/coordinate$quality/${sample}-$region.Combined.tab | awk '$2 == "T" || $2 == "t"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-Ribo.$region.$nuc.bed

		elif [[ $nuc == "Combined" ]]; then
			echo $A_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-RiboCombined.$region.A.txt
			echo $C_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-RiboCombined.$region.C.txt
			echo $G_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-RiboCombined.$region.G.txt
			echo $U_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-RiboCombined.$region.U.txt
				
			#Create BED file for all ribonucleotides
			cp $repository/results/$sample/coordinate$quality/${sample}-$region.Combined.tab $output/${sample}-Ribo.$region.$nuc.bed
		fi
		
#############################################################################################################################
			
		#Obtain coordinates/sequences of dNMPs +/- 100 bp from rNMPs
		
		#Continue only for BED files > 0
		if [[ -s $output/${sample}-Ribo.$region.$nuc.bed ]]; then
			
			#Obtain coordinates of flanking sequences and remove coordinates where start = end
			bedtools flank -i $output/${sample}-Ribo.$region.$nuc.bed -s -g $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes -l 100 -r 0 | awk '$2 != $3' > $output/${sample}-Upstream.$region.$nuc.bed
			bedtools flank -i $output/${sample}-Ribo.$region.$nuc.bed -s -g $(dirname $fasta)/$(basename $fasta .fa).chrom.sizes -l 0 -r 100 | awk '$2 != $3' > $output/${sample}-Downstream.$region.$nuc.bed
	
			#Obtain nucleotides flanking rNMPs (reverse order of upstream) and insert tabs bases for easier parsing
			bedtools getfasta -s -fi $repository/$(basename $fasta .fa)-$region.fa -bed $output/${sample}-Upstream.$region.$nuc.bed | grep -v '>' | rev | sed 's/.../& /2g;s/./& /g' > $output/${sample}-Upstream.$region.$nuc.tab
			bedtools getfasta -s -fi $repository/$(basename $fasta .fa)-$region.fa -bed $output/${sample}-Downstream.$region.$nuc.bed | grep -v '>' | sed 's/.../& /2g;s/./& /g' > $output/${sample}-Downstream.$region.$nuc.tab
							
			#Calculate frequencies of dNMPs +/- 100 base pairs from rNMPs
			for direction in "Upstream" "Downstream"; do
		
				for i in {1..100}; do
		
					#Calculate count of each dNMP
					A_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$region.$nuc.tab | grep -Eo 'A|a' | wc -l)
					C_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$region.$nuc.tab | grep -Eo 'C|c' | wc -l)
					G_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$region.$nuc.tab | grep -Eo 'G|g' | wc -l)
					T_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$region.$nuc.tab | grep -Eo 'T|t' | wc -l)

					#Calculate total number of dNMPs
					FlankTotal=$(($A_Flank + $C_Flank + $G_Flank + $T_Flank))

					#Calculate raw frequency of each dNMP
					if [[ $FlankTotal == 0 ]]; then
						
						A_FlankFreq='NA'
						C_FlankFreq='NA'
						G_FlankFreq='NA'
						T_FlankFreq='NA'
				
					else
						
						A_FlankFreq=$(echo "($A_Flank/$FlankTotal)" | bc -l)
						C_FlankFreq=$(echo "($C_Flank/$FlankTotal)" | bc -l)
						G_FlankFreq=$(echo "($G_Flank/$FlankTotal)" | bc -l)
						T_FlankFreq=$(echo "($T_Flank/$FlankTotal)" | bc -l)
					fi
				
					#Save dNMPs frequencies to .txt files
					if [[ $A_FlankFreq == 'NA' ]]; then
						echo $A_FlankFreq >> $output/${sample}-$direction.$region.$nuc.A.txt
					else
						echo $A_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.$region.$nuc.A.txt
					fi
				
					if [[ $C_FlankFreq != 'NA' ]]; then
						echo $C_FlankFreq >> $output/${sample}-$direction.$region.$nuc.C.txt
					else
						echo $C_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.$region.$nuc.C.txt
					fi
				
					if [[ $G_FlankFreq != 'NA' ]]; then
						echo $G_FlankFreq >> $output/${sample}-$direction.$region.$nuc.G.txt
					else
						echo $G_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.$region.$nuc.G.txt	
					fi
				
					if [[ $T_FlankFreq != 'NA' ]]; then
						echo $T_FlankFreq >> $output/${sample}-$direction.$region.$nuc.T.txt
					else
						echo $T_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.$region.$nuc.T.txt
					fi

				done
			done

#############################################################################################################################
			
			#Create and save dataset file containing nucleotide frequencies
			
			#Add nucleotides to header line
			echo -e "\tA\tC\tG\tU/T" > $output/${sample}-$region.$nuc.raw.tab
			echo -e "\tA\tC\tG\tU/T" > $output/${sample}-$region.$nuc.normalized.tab
				
			#Background frequencies of dNMPs
			A_BkgFreq=$(cut -f2 $repository/$(basename $fasta .fa)-$region.txt | head -1)
			C_BkgFreq=$(cut -f2 $repository/$(basename $fasta .fa)-$region.txt | head -2 | tail -1)
			G_BkgFreq=$(cut -f2 $repository/$(basename $fasta .fa)-$region.txt | head -3 | tail -1)
			T_BkgFreq=$(cut -f2 $repository/$(basename $fasta .fa)-$region.txt | head -4 | tail -1)
			
			#Save rNMP and flanking frequencies to separate variables
			A_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.$region.$nuc.A.txt | tac) <(cat $output/${sample}-Ribo$nuc.$region.A.txt) <(cat $output/${sample}-Downstream.$region.$nuc.A.txt)))
			C_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.$region.$nuc.C.txt | tac) <(cat $output/${sample}-Ribo$nuc.$region.C.txt) <(cat $output/${sample}-Downstream.$region.$nuc.C.txt)))
			G_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.$region.$nuc.G.txt | tac) <(cat $output/${sample}-Ribo$nuc.$region.G.txt) <(cat $output/${sample}-Downstream.$region.$nuc.G.txt)))
			T_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.$region.$nuc.T.txt | tac) <(cat $output/${sample}-Ribo$nuc.$region.U.txt) <(cat $output/${sample}-Downstream.$region.$nuc.T.txt)))

			#Add positions and frequencies of nucleotides in correct order to create dataset (un-normalized to reference genome)
			paste <(echo "$(seq -100 1 100)") <(for a in $A_allFreqs; do if [[ $a != "NA" ]]; then echo $a | bc -l | xargs printf "%.*f\n" 5; else echo $a; fi; done) \
							  <(for c in $C_allFreqs; do if [[ $c != "NA" ]]; then echo $c | bc -l | xargs printf "%.*f\n" 5; else echo $c; fi; done) \
							  <(for g in $G_allFreqs; do if [[ $g != "NA" ]]; then echo $g | bc -l | xargs printf "%.*f\n" 5; else echo $g; fi; done) \
							  <(for t in $T_allFreqs; do if [[ $t != "NA" ]]; then echo $t | bc -l | xargs printf "%.*f\n" 5; else echo $t; fi; done) \
							  >> $output/${sample}-$region.$nuc.raw.tab
								  
			#Add positions and frequencies of nucleotides in correct order to create dataset (normalized to reference genome)
			paste <(echo "$(seq -100 1 100)") <(for a in $A_allFreqs; do if [[ $a != "NA" ]]; then echo $a/$A_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $a; fi; done) \
							  <(for c in $C_allFreqs; do if [[ $c != "NA" ]]; then echo $c/$C_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $c; fi; done) \
							  <(for g in $G_allFreqs; do if [[ $g != "NA" ]]; then echo $g/$G_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $g; fi; done) \
							  <(for t in $T_allFreqs; do if [[ $t != "NA" ]]; then echo $t/$T_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $t; fi; done) \
							  >> $output/${sample}-$region.$nuc.normalized.tab

#############################################################################################################################
			#Print status
			echo "Status: Sequence Module for $sample ($nuc,$region) is complete"
				
		fi
	done
	
	#Remove temporary files
	rm $output/*.txt $output/*{Upstream,Downstream}*.{bed,tab}

done

rm $output/*.bed
