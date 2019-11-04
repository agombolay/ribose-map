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
	paste <(echo -e "A") <(echo "$A_BkgFreq" | xargs printf "%.*f\n" 5) > $(dirname $fasta)/$(basename $fasta .fa)-$region.txt
	paste <(echo -e "C") <(echo "$C_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.txt
	paste <(echo -e "G") <(echo "$G_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.txt
	paste <(echo -e "U") <(echo "$T_BkgFreq" | xargs printf "%.*f\n" 5) >> $(dirname $fasta)/$(basename $fasta .fa)-$region.txt

done

#############################################################################################################################

#Calculate frequencies of rNMPs
for file in $output/${sample}-*.bed; do

	if [ -s $file ]; then
	
		temp=$(echo $file | awk -F '[-]' '{print $2 $3 $4}')
		region=$(basename $temp .bed)
		
		#Extract rNMP nucleotides from FASTA
		bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -bed $file | grep -v '>' > $output/{$sample}-$region.ribos.txt
			
		#Calculate counts of rNMPs
		A_Ribo=$(awk '$1 == "A" || $1 == "a"' $output/{$sample}-$region.ribos.txt | wc -l)
		C_Ribo=$(awk '$1 == "C" || $1 == "c"' $output/{$sample}-$region.ribos.txt | wc -l)
		G_Ribo=$(awk '$1 == "G" || $1 == "g"' $output/{$sample}-$region.ribos.txt | wc -l)
		U_Ribo=$(awk '$1 == "T" || $1 == "t"' $output/{$sample}-$region.ribos.txt | wc -l)
	
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
				echo $A_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.G_Ribo.txt
				echo 'NA' > $output/${sample}-$region.U_Ribo.txt
				
				#Create BED file for only A ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "A" || $2 == "a"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.A.bed

			elif [[ $nuc == "C" ]]; then
				echo $C_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.G_Ribo.txt
				echo 'NA' > $output/${sample}-$region.U_Ribo.txt
				
				#Create BED file for only C ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "C" || $2 == "c"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.C.bed
				
			elif [[ $nuc == "G" ]]; then
				echo $G_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.G_Ribo.txt
				echo 'NA' > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.U_Ribo.txt
				
				#Create BED file for only G ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "G" || $2 == "g"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.G.bed

			elif [[ $nuc == "U" ]]; then
				echo $U_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.U_Ribo.txt
				echo 'NA' > $output/${sample}-$region.A_Ribo.txt
				echo 'NA' > $output/${sample}-$region.C_Ribo.txt
				echo 'NA' > $output/${sample}-$region.G_Ribo.txt
				
				#Create BED file for only U ribonucleotide
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)-$region.fa -tab -bed $file | awk '$2 == "T" || $2 == "t"' | cut -f1 | sed 's/\:/\t/' | sed 's/\-/\t/' | sed 's/(/\t.\t.\t/;s/)//' > $output/${sample}-$region.U.bed

			elif [[ $nuc == "Combined" ]]; then
				echo $A_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.A_Ribo.txt
				echo $C_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.C_Ribo.txt
				echo $G_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.G_Ribo.txt
				echo $U_RiboFreq | xargs printf "%.*f\n" 5 > $output/${sample}-$region.U_Ribo.txt
				
				#Create BED file for all ribonucleotides
				cp $file $output/${sample}-$region.Combined.bed
			fi
		
#############################################################################################################################
			
			#Obtain coordinates/sequences of dNMPs +/- 100 bp from rNMPs
		
			#Continue only for BED files > 0
			if [[ -s $output/${sample}-$region.$nuc.bed ]]; then
			
				#Obtain coordinates of flanking sequences and remove coordinates where start = end
				bedtools flank -i $output/${sample}-$region.$nuc.bed -s -g $(dirname $fasta)/$(basename $fasta .fa).bed -l 100 -r 0 | awk '$2 != $3' > $output/${sample}-Upstream.$nuc.bed
				bedtools flank -i $output/${sample}-$region.$nuc.bed -s -g $(dirname $fasta)/$(basename $fasta .fa).bed -l 0 -r 100 | awk '$2 != $3' > $output/${sample}-Downstream.$nuc.bed
	
				#Obtain nucleotides flanking rNMPs (reverse order of upstream) and insert tabs bases for easier parsing
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)_$region.fa -bed $output/Upstream.bed | grep -v '>' | rev | sed 's/.../& /2g;s/./& /g' > $output/${sample}-Upstream.$nuc.tab
				bedtools getfasta -s -fi $(dirname $fasta)/$(basename $fasta .fa)_$region.fa -bed $output/Downstream.bed | grep -v '>' | sed 's/.../& /2g;s/./& /g' > $output/${sample}-Downstream.$nuc.tab
							
				#Calculate frequencies of dNMPs +/- 100 base pairs from rNMPs
				for direction in "Upstream" "Downstream"; do
		
					for i in {1..100}; do
		
						#Calculate count of each dNMP
						A_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$nuc.tab | grep -Eo 'A|a' | wc -l)
						C_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$nuc.tab | grep -Eo 'C|c' | wc -l)
						G_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$nuc.tab | grep -Eo 'G|g' | wc -l)
						T_Flank=$(awk -v field=$i '{ print $field }' $output/${sample}-$direction.$nuc.tab | grep -Eo 'T|t' | wc -l)

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
							echo $A_FlankFreq >> $output/${sample}-$direction.A.txt
						else
							echo $A_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.A.txt
						fi
				
						if [[ $C_FlankFreq != 'NA' ]]; then
							echo $C_FlankFreq >> $output/${sample}-$direction.C.txt
						else
							echo $C_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.C.txt
						fi
				
						if [[ $G_FlankFreq != 'NA' ]]; then
							echo $G_FlankFreq >> $output/${sample}-$direction.G.txt
						else
							echo $G_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.G.txt	
						fi
				
						if [[ $T_FlankFreq != 'NA' ]]; then
							echo $T_FlankFreq >> $output/${sample}-$direction.T.txt
						else
							echo $T_FlankFreq | xargs printf "%.*f\n" 5 >> $output/${sample}-$direction.T.txt
						fi

					done
				done

#############################################################################################################################
				#Create and save dataset file containing nucleotide frequencies
			
				#Add nucleotides to header line
				echo -e "\tA\tC\tG\tU/T" > $output/${sample}-$region.$nuc.raw.tab
				echo -e "\tA\tC\tG\tU/T" > $output/${sample}-$region.$nuc.normalized.tab
				
				#Background frequencies of dNMPs
				A_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -1)
				C_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -2 | tail -1)
				G_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -3 | tail -1)
				T_BkgFreq=$(cut -f2 $(dirname $fasta)/$(basename $fasta .fa)-$region.txt | head -4 | tail -1)
			
				#Save rNMP and flanking frequencies to separate variables
				A_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.A.txt | tac) <(cat $output/${sample}-$region.A_Ribo.txt) <(cat $output/${sample}-Downstream.A.txt)))
				C_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.C.txt | tac) <(cat $output/${sample}-$region.C_Ribo.txt) <(cat $output/${sample}-Downstream.C.txt)))
				G_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.G.txt | tac) <(cat $output/${sample}-$region.G_Ribo.txt) <(cat $output/${sample}-Downstream.G.txt)))
				T_allFreqs=$(paste <(cat <(cat $output/${sample}-Upstream.T.txt | tac) <(cat $output/${sample}-$region.U_Ribo.txt) <(cat $output/${sample}-Downstream.T.txt)))

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
	fi
	
	#Remove temporary files
	#rm $output/*.txt $output/*{Upstream,Downstream}.{bed,tab}

done
