background="/projects/home/agombolay3/data/repository/Ribose-seq-Project/ribose-seq/reference/sacCer2.fa"
nucleotides="/projects/home/agombolay3/data/repository/Ribose-seq-Project/ribose-seq/results/sacCer2/FS15.trimmed/Nucleotide-Frequencies/Nucleotides/Columns/nuclear/upstream/*.txt"

for file in $background;
do
  	A=$(grep -v '>' $file | grep -o 'A' - | wc -l)
        C=$(grep -v '>' $file | grep -o 'C' - | wc -l)
        G=$(grep -v '>' $file | grep -o 'G' - | wc -l)
        T=$(grep -v '>' $file | grep -o 'T' - | wc -l)

        total=$(($A+$C+$G+$T))

        A_background_frequency=$(bc <<< "scale = 4; `expr $A/$total`")
        C_background_frequency=$(bc <<< "scale = 4; `expr $C/$total`")
        G_background_frequency=$(bc <<< "scale = 4; `expr $G/$total`")
        T_background_frequency=$(bc <<< "scale = 4; `expr $T/$total`")
	
	echo "A Background Frequency: $A_background_frequency" >> background.txt
	echo "C Background Frequency: $C_background_frequency" >> background.txt
	echo "G Background Frequency: $G_background_frequency" >> background.txt
	echo "T Background Frequency: $T_background_frequency" >> background.txt

done

for file in $nucleotides;
do
	A=$(grep -v '>' $file | grep -o 'A' - | wc -l)
	C=$(grep -v '>' $file | grep -o 'C' - | wc -l)
	G=$(grep -v '>' $file | grep -o 'G' - | wc -l)
	T=$(grep -v '>' $file | grep -o 'T' - | wc -l)

	total=$(($A+$C+$G+$T))
	
	A_frequency=$(bc <<< "scale = 4; `expr $A/$total`")
	C_frequency=$(bc <<< "scale = 4; `expr $C/$total`")
	G_frequency=$(bc <<< "scale = 4; `expr $G/$total`")
	T_frequency=$(bc <<< "scale = 4; `expr $T/$total`")

	A_normalized_frequency=$(bc <<< "scale = 4; `expr $A_frequency/$A_background_frequency`")
        C_normalized_frequency=$(bc <<< "scale = 4; `expr $C_frequency/$C_background_frequency`")
        G_normalized_frequency=$(bc <<< "scale = 4; `expr $G_frequency/$G_background_frequency`")
        T_normalized_frequency=$(bc <<< "scale = 4; `expr $T_frequency/$T_background_frequency`")

	echo $A_normalized_frequency >> A_normalized_frequencies.txt
	echo $C_normalized_frequency >> C_normalized_frequencies.txt
	echo $G_normalized_frequency >> G_normalized_frequencies.txt
	echo $T_normalized_frequency >> T_normalized_frequencies.txt

	paste A_normalized_frequencies.txt C_normalized_frequencies.txt G_normalized_frequencies.txt T_normalized_frequencies.txt > normalized_frequencies.txt

done
