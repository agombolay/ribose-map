#Count number of reads with correct barcode
samtools view -h sample.bam -o sample.sam
grep -e 'UMI_....barcode....' -e '@HG' -e 'SQ' -e 'PG' sample.sam > sample-filtered.sam
samtools view sample-filtered.sam -b -S | samtools sort -o sample-filtered.bam
samtools index sample-filtered.bam
