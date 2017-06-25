<h1><p align="center">S. pombe (positive control for ribose-seq)</p></h1>

<img src="https://github.com/agombolay/Ribose-seq-Project/blob/master/Additional-Information/pombe.png?raw=true" width="300px" height="100px" />

**Homology Domain H1 of mat1 locus (Diribo Imprint: chrII:2115362-63 on forward strand)**:
TAA**TT**TTTTTGTAATATAAATGTATAGTCTTTCTCCTTTGTTTTCTCTCGTTCGTTTCCATGT

## Commands to determine coordinates of imprint
```
bowtie2 -x pombe -c -U TAATTTTTTTGTAATATAAATGTATAGTCTTTCTCCTTTGTTTTCTCTCGTTCGTTTCCATGT -S pombe.sam
samtools view pombe.sam -b | samtools sort - -o pombe.bam
samtools index pombe.bam
```
```
samtools tview -p II:2115359 pombe.bam /projects/home/agombolay3/data/repository/Ribose-seq-Project/ribose-seq/reference/pombe.fa
```

**Reference**:
Vengrova and Dalgaard, Genes and Development 2004

**FASTA Sequence of S. pombe genome**:  
http://fungi.ensembl.org/info/website/ftp/index.html
