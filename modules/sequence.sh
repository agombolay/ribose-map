#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#Calculate frequencies of rNMP nucleotides and those of flanking nucleotides

#############################################################################################################################
#Load config file
. "$1"

#Create output directory and remove any old files
output=$repository/results/$sample/sequence$quality
rm -r $output; mkdir -p $output

#############################################################################################################################

#Create .fai file for reference
samtools faidx $fasta

#Create .bed file for reference
cut -f 1,2 $fasta.fai > $(dirname $fasta)/$(basename $fasta .fa).bed

#############################################################################################################################

if [[ $mito ]]; then
	
	if [[ $other ]]; then
	
		other_new=$(echo $other | sed 's/ /|/g')
		
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito - | grep -Ewv $other_new -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa

		#Extract only unique coordinates and then subset based on region
		awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab | grep -wv $mito - | grep -Ewv $other_new - > $output/${sample}-nucleus.bed
		
		#Mito
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-mito.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-mito.fa

		#Extract only unique coordinates and then subset based on region
		awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab | grep -w $mito - > $output/${sample}-mito.bed
		
		#Other
		for region in $other; do
			
			#Subset FASTA file based on region
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			
			#Extract only unique coordinates and then subset based on region
			awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab | grep -w $region - > $output/${sample}-$region.bed
		
		done

	else
	
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -wv $mito -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		
		#Extract only unique coordinates and then subset based on region
		awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab | grep -wv $mito - > $output/${sample}-nucleus.bed

		#Mito
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $mito -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-mito.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-mito.fa
		
		#Extract only unique coordinates and then subset based on region
		awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab | grep -w $mito - > $output/${sample}-mito.bed

	fi
	
elif [[ ! $mito ]]; then
		
	if [[ $other ]]; then
		
		other_new=$(echo $other | sed 's/ /|/g')
		
		#Nucleus
		#Subset FASTA file based on region
		chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -Ewv $other_new -)
		samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		
		#Extract only unique coordinates and then subset based on region
		awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab | grep -Ewv $other_new - > $output/${sample}-nucleus.bed

		#Other
		for region in $other; do
			
			#Subset FASTA file based on region
			chr=$(awk '{print $1}' $(dirname $fasta)/$(basename $fasta .fa).bed | grep -w $region -)
			samtools faidx $fasta $chr > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
			
			#Extract only unique coordinates and then subset based on region
			awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab | grep -w $region - > $output/${sample}-$region.bed

		done

	else
	
		#Nucleus
		#Subset FASTA file based on region
		cat $fasta > $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		samtools faidx $(dirname $fasta)/$(basename $fasta .fa)-nucleus.fa
		
		#Extract only unique coordinates and then subset based on region
		awk -v "OFS=\t" '{print $1, $2, $3, ".", ".", $4}' $repository/results/$sample/coordinate$quality/$sample.counts.tab > $output/${sample}-nucleus.bed

	fi

fi

#############################################################################################################################

#Calculate frequencies of reference genome
for file in $(dirname $fasta)/$(basename $fasta .fa)-*.fa; do
			
	temp=$(echo $file | awk -F '[-]' '{print $2 $3 $4}')
	region=$(basename $temp .fa)
			
	#Nucleotide Frequencies of Reference Genome
			
	#Calculate counts of each nucleotide
	A_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)-$region.fa | grep -Eo 'A|a' - | wc -l)
	C_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)-$region.fa | grep -Eo 'C|c' - | wc -l)
	G_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)-$region.fa | grep -Eo 'G|g' - | wc -l)
	T_Bkg=$(grep -v '>' $(dirname $fasta)/$(basename $fasta .fa)-$region.fa | grep -Eo 'T|t' - | wc -l)
	
	#Calculate total number of nucleotides
	BkgTotal=$(($A_Bkg + $C_Bkg + $G_Bkg + $T_Bkg))
		
	#Calculate frequencies of each nucleotide
	A_BkgFreq=$(echo "($A_Bkg + $T_Bkg)/($BkgTotal*2)" | bc -l)
	C_BkgFreq=$(echo "($C_Bkg + $G_Bkg)/($BkgTotal*2)" | bc -l)
	G_BkgFreq=$(echo "($G_Bkg + $C_Bkg)/($BkgTotal*2)" | bc -l)
	T_BkgFreq=$(echo "($T_Bkg + $A_Bkg)/($BkgTotal*2)" | bc -l)
			
	#Save nucleotide frequencies to .txt file
	paste <(echo -e "A") <(echo "$A_BkgFreq" | xargs printf "%.*f\n" 5) > $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
	paste <(echo -e "C") <(echo "$C_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
	paste <(echo -e "G") <(echo "$G_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.fa
	paste <(echo -e "U") <(echo "$T_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.fa

done

#############################################################################################################################

#Calculate frequencies of rNMPs
for file in $output/${sample}-*.bed; do

	if [ -s $output/${sample}-*.bed ]; then
	
		temp=$(echo $file | awk -F '[-]' '{print $2 $3 $4}')
		region=$(basename $temp .bed)
		
		for nuc in "A" "C" "G" "T" "Combined"; do
			
			#Extract rNMP nucleotides from FASTA
			bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -bed $file | grep -v '>' > $output/{$sample}-$region.ribos.txt
			
			#Calculate counts of rNMPs
			A_Ribo=$(awk '$1 == "A" || $1 == "a"' $output/{$sample}-$region.ribos.txt | wc -l)
			C_Ribo=$(awk '$1 == "C" || $1 == "c"' $output/{$sample}-$region.ribos.txt | wc -l)
			G_Ribo=$(awk '$1 == "G" || $1 == "g"' $output/{$sample}-$region.ribos.txt | wc -l)
			U_Ribo=$(awk '$1 == "T" || $1 == "t"' $output/{$sample}-$region.ribos.txt | wc -l)
	
			#Calculate total number of rNMPs
			RiboTotal=$(($A_Ribo + $C_Ribo + $G_Ribo + $U_Ribo))
	
			#Calculate normalized frequency of each rNMP
			A_RiboFreq=$(echo "($A_Ribo/$RiboTotal)/$A_BkgFreq" | bc -l)
			C_RiboFreq=$(echo "($C_Ribo/$RiboTotal)/$C_BkgFreq" | bc -l)
			G_RiboFreq=$(echo "($G_Ribo/$RiboTotal)/$G_BkgFreq" | bc -l)
			U_RiboFreq=$(echo "($U_Ribo/$RiboTotal)/$T_BkgFreq" | bc -l)
				
			#Calculate un-normalized frequency of each rNMP
			#A_RiboFreq=$(echo "($A_Ribo/$RiboTotal)" | bc -l)
			#C_RiboFreq=$(echo "($C_Ribo/$RiboTotal)" | bc -l)
			#G_RiboFreq=$(echo "($G_Ribo/$RiboTotal)" | bc -l)
			#U_RiboFreq=$(echo "($U_Ribo/$RiboTotal)" | bc -l)
			
			#Save frequencies of rNMPs to TXT files
			if [[ $nuc == "A" ]]; then
				echo $A_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.G_Ribo.txt
				echo 'NA' > $output/${sample}-$region.U_Ribo.txt
				
			elif [[ $nuc == "C" ]]; then
				echo $C_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.G_Ribo.txt
				echo 'NA' > $output/${sample}-$region.U_Ribo.txt
				
			elif [[ $nuc == "G" ]]; then
				echo $G_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.G_Ribo.txt
				echo 'NA' > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.U_Ribo.txt
				
			elif [[ $nuc == "T" ]]; then
				echo $U_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.U_Ribo.txt
				echo 'NA' > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.G_Ribo.txt
				
			elif [[ $nuc == "Combined" ]]; then
				echo $A_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.A_Ribo.txt
				echo $C_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.C_Ribo.txt
				echo $G_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.G_Ribo.txt
				echo $U_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.U_Ribo.txt
			fi

			#Combine rNMP frequencies into one file
			Ribo=$(paste $output/${sample}-$region.{A,C,G,U}_Ribo.txt)

#############################################################################################################################
			
			#Obtain coordinates/sequences of dNMPs +/- 100 bp from rNMPs

			#Create 5 BED files, one for each nucleotide and one combined
			if [[ $nuc == "A" ]]; then
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "A" || $2 == "a"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.$nuc.bed
			
			elif [[ $nuc == "C" ]]; then
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "C" || $2 == "c"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.$nuc.bed
			
			elif [[ $nuc == "G" ]]; then
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "G" || $2 == "g"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.$nuc.bed
			
			elif [[ $nuc == "T" ]]; then
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "T" || $2 == "t"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.$nuc.bed
			
			elif [[ $nuc == "Combined" ]]; then
				cp $file $output/${sample}-$region.$nuc.bed 
			fi
		
			#Continue only for BED files > 0
			if [[ -s $output/${sample}-$region.$nuc.bed ]]; then
			
				#Obtain coordinates of flanking sequences and remove coordinates where start = end
				bedtools flank -i $output/${sample}-$region.$nuc.bed -s -g $(dirname $fasta)/$(basename $fasta .fa).bed -l 100 -r 0 | awk '$2 != $3' > $output/Up.bed
				bedtools flank -i $output/${sample}-$region.$nuc.bed -s -g $(dirname $fasta)/$(basename $fasta .fa).bed -l 0 -r 100 | awk '$2 != $3' > $output/Down.bed
	
				#Obtain nucleotides flanking rNMPs (reverse order of up) and insert tabs bases for easier parsing
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)_$region.fa -bed $output/Down.bed | grep -v '>' | sed 's/.../& /2g;s/./& /g' > $output/Down.tab
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)_$region.fa -bed $output/Up.bed | grep -v '>' | rev | sed 's/.../& /2g;s/./& /g' > $output/Up.tab
				
#############################################################################################################################
				#Calculate frequencies of dNMPs +/- 100 base pairs from rNMPs

				for dir in "Up" "Down"; do
		
					for i in {1..100}; do
		
						#Calculate count of each dNMP
						A_Flank=$(awk -v field=$i '{ print $field }' $output/$dir.tab | grep -Eo 'A|a' | wc -l)
						C_Flank=$(awk -v field=$i '{ print $field }' $output/$dir.tab | grep -Eo 'C|c' | wc -l)
						G_Flank=$(awk -v field=$i '{ print $field }' $output/$dir.tab | grep -Eo 'G|g' | wc -l)
						T_Flank=$(awk -v field=$i '{ print $field }' $output/$dir.tab | grep -Eo 'T|t' | wc -l)

						#Calculate total number of dNMPs
						FlankTotal=$(($A_Flank + $C_Flank + $G_Flank + $T_Flank))

						#Calculate normalized frequencies of dNMPs
						if [[ $FlankTotal != 0 ]]; then
						
							A_FlankFreq=$(echo "($A_Flank/$FlankTotal)/$A_BkgFreq" | bc -l)
							C_FlankFreq=$(echo "($C_Flank/$FlankTotal)/$C_BkgFreq" | bc -l)
							G_FlankFreq=$(echo "($G_Flank/$FlankTotal)/$G_BkgFreq" | bc -l)
							T_FlankFreq=$(echo "($T_Flank/$FlankTotal)/$T_BkgFreq" | bc -l)
								
							#A_FlankFreq=$(echo "($A_Flank/$FlankTotal)" | bc -l)
							#C_FlankFreq=$(echo "($C_Flank/$FlankTotal)" | bc -l)
							#G_FlankFreq=$(echo "($G_Flank/$FlankTotal)" | bc -l)
							#T_FlankFreq=$(echo "($T_Flank/$FlankTotal)" | bc -l)
				
						elif [[ $FlankTotal == 0 ]]; then
							
							A_FlankFreq='NA'
							C_FlankFreq='NA'
							G_FlankFreq='NA'
							T_FlankFreq='NA'
						fi
				
						#Save normalized dNMPs frequencies to TXT files
						if [[ $A_FlankFreq != 'NA' ]]; then
							echo $A_FlankFreq | xargs printf "%.*f\n" 5 >> $output/A_$dir.txt
						
						elif [[ $A_FlankFreq == 'NA' ]]; then
							echo $A_FlankFreq >> $output/A_$dir.txt
						fi
				
						if [[ $C_FlankFreq != 'NA' ]]; then
							echo $C_FlankFreq | xargs printf "%.*f\n" 5 >> $output/C_$dir.txt
						
						elif [[ $C_FlankFreq == 'NA' ]]; then
							echo $C_FlankFreq >> $output/C_$dir.txt
						fi
				
						if [[ $G_FlankFreq != 'NA' ]]; then
							echo $G_FlankFreq | xargs printf "%.*f\n" 5 >> $output/G_$dir.txt
						
						elif [[ $G_FlankFreq == 'NA' ]]; then
							echo $G_FlankFreq >> $output/G_$dir.txt
						fi
				
						if [[ $T_FlankFreq != 'NA' ]]; then
							echo $T_FlankFreq | xargs printf "%.*f\n" 5 >> $output/T_$dir.txt
						
						elif [[ $T_FlankFreq == 'NA' ]]; then
							echo $T_FlankFreq >> $output/T_$dir.txt
						fi
		
						#Combine dNMP frequencies into one file per location
						if [[ $dir == "Up" ]]; then
							#Print upstream frequencies in reverse order
							Up=$(paste $output/{A,C,G,T}_Up.txt | tac -)
						
						elif [[ $dir == "Down" ]]; then
							Down=$(paste $output/{A,C,G,T}_Down.txt)
						fi

					done
				done

#############################################################################################################################
					#Create and save dataset file containing nucleotide frequencies
			
					#Add nucleotides to header line
					#echo -e "\tA\tC\tG\tU/T" > $output/$sample.$nuc.$region.raw.tab
					echo -e "\tA\tC\tG\tU/T" > $output/$sample.$nuc.$region.normalized.tab
			
					A=$(paste <(cat <(cat $output/A_Up.txt | tac) <(cat $output/A_Ribo.txt) <(cat $output/A_Down.txt)))
					C=$(paste <(cat <(cat $output/C_Up.txt | tac) <(cat $output/C_Ribo.txt) <(cat $output/C_Down.txt)))
					G=$(paste <(cat <(cat $output/G_Up.txt | tac) <(cat $output/G_Ribo.txt) <(cat $output/G_Down.txt)))
					T=$(paste <(cat <(cat $output/T_Up.txt | tac) <(cat $output/U_Ribo.txt) <(cat $output/T_Down.txt)))

					#Add positions and frequencies of nucleotides in correct order to create dataset (not normalized to anything)
					#paste <(echo "$(seq -100 1 100)") <(cat <(echo "$Up") <(echo "$Ribo") <(echo "$Down")) >> $output/$sample.$nuc.$region.raw.tab
					
					#Add positions and frequencies of nucleotides in correct order to create dataset (normalized to reference genome)
					paste <(echo "$(seq -100 1 100)") <(for i in $A; do if [[ $i != "NA" ]]; then echo $i/$A_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $i; fi; done) \
									  <(for j in $C; do if [[ $j != "NA" ]]; then echo $j/$C_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $j; fi; done) \
									  <(for k in $G; do if [[ $k != "NA" ]]; then echo $k/$G_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $k; fi; done) \
									  <(for l in $T; do if [[ $l != "NA" ]]; then echo $l/$T_BkgFreq | bc -l | xargs printf "%.*f\n" 5; else echo $l; fi; done) \
									  >> $output/$sample.$nuc.$region.normalized.tab

#############################################################################################################################
					#Print status
					echo "Status: Sequence Module for $sample ($nuc,$region) is complete"
					
				fi
			fi
		
			#Remove temporary files
			rm -f $output/*.{txt,fa,fa.fai} $output/{Up,Down}.bed $output/{Up,Down}.tab
	
	done
done
