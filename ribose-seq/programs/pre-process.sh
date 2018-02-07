#Extract UMI from 5' ends of reads
if [[ $umi ]] && [[ ! $read2 ]]; then
	umi_tools extract -v 0 -I $fastq1 --bc-pattern=$UMI -S processed.fq.gz
	
elif [[ $umi ]] && [[ $read2 ]]; then
	umi_tools extract -v 0 -I $fastq1 --bc-pattern=$UMI --read2-in=$fastq2 \
	--stdout=processed1.fq.gz --read2-out=processed2.fq.gz
fi

#############################################################################################################################
if [[ $barcode ]]; then
	#Filter FASTQ file based on barcode sequence
	grep --no-group-separator -B1 -A2 ^$barcode $output/UMI.fq > $output/filtered.fq
fi

#############################################################################################################################
if [[ ! $read2 ]]; then
	if [[ ! $adapter ]]; then
		if [[ ! $barcode ]]; then
			trim_galore --gzip --length $min $output/file1.fq -o $output
		elif [[ $barcode ]]; then
			grep --no-group-separator -B1 -A2 ^$barcode $output/file1.fq \
			| trim_galore --gzip --length $min --clip_R1 3 - -o $output
			
	elif [[ $adapter ]]; then
		if [[ ! $barcode ]]; then
			trim_galore --gzip --length $min -a $adapter $output/file1.fq -o $output
		elif [[ $barcode ]]; then
			grep --no-group-separator -B1 -A2 ^$barcode $output/file1.fq | trim_galore \
			--gzip--length $min --clip_R1 3 -a $adapter - -o $output

elif [[ $read2 ]]; then
	
	if [[ ! $adapter ]]; then
		if [[ ! $barcode ]]; then
			trim_galore --gzip --paired --length $min $output/processed1.fq -o $output
		elif [[ $barcode ]]; then
			grep --no-group-separator -B1 -A2 ^$barcode $output/processed1.fq > $output/processed2.fq
			trim_galore --gzip --paired --length $minimum --clip_R1 3 $output/filtered.fq -o $output
			
	elif [[ $adapter ]]; then
		if [[ ! $barcode ]]; then
			trim_galore --gzip --paired --length $minimum -a $adapter $output/processed1.fq -o $output
		elif [[ $barcode ]]; then
			grep --no-group-separator -B1 -A2 ^$barcode $output/processed1.fq > $output/processed2.fq
			trim_galore --gzip --paired --length $minimum --clip_R1 3 -a $adapter $output/filtered.fq -o $output

fi
