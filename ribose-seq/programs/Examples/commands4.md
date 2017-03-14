#Count number of reads with correct barcode
```
samtools view -h sample.bam -o sample.sam
grep -e 'UMI_....barcode....' -e '@HG' -e '@SQ' -e '@PG' sample.sam > sample-filtered.sam
samtools view sample-filtered.sam -b -S | samtools sort -o sample-filtered.bam
samtools index sample-filtered.bam

rm sample.sam sample-filtered.sam
```

Results of these commands should match!
```
samtools view -c -b sample-filtered.bam
samtools bam2fq sample.bam | grep 'UMI_....barcode....' - | wc -l
```
