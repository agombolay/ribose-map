bedtools bamtofastq -i sample.bam -fq sample.fastq
grep 'UMI_....TGA....' sample.fastq | wc -l 
