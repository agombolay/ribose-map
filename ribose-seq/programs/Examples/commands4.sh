#Count number of reads with correct barcode
samtools view -h sample.bam -o sample.sam
grep 'UMI_....barcode....' sample.fastq | wc -l
samtools view sample.sam -b -S | samtools sort -o sample.bam
samtools index sample.bam
