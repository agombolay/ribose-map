<h1><p align="center">Count number of reads</p></h1>
* [Count number of reads in a FASTA file] (http://thegenomefactory.blogspot.com/2011/09/counting-sequences-with-unix-tools.html)

* [Count number of reads in BAM file] (http://crazyhottommy.blogspot.com/2013/06/count-how-many-mapped-reads-in-bam-file.html)

**Mapped reads**:
```
samtools view -c -F 4
```

**Unmapped reads**:
```
samtools view -c -f 4
```

<h1><p align="center">Check quality of reads</p></h1>
* [FastQC] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

<h1><p align="center">Trim reads based on quality and adapters</p></h1>
Importance of removing adapters: Increase mapping percentage and decrease incorrect mappings

* [Trimmomatic] (http://www.usadellab.org/cms/?page=trimmomatic)
  * [Manual for Trimmomatic] (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

<h1><p align="center">View alignment data</p></h1>

###[SAMtools tview] (http://samtools.sourceforge.net/tview.shtml)

```
samtools tview <BAM> --reference <Fasta>
```

<h1><p align="center">Convert file formats</p></h1>

###[BAM to FASTA:] (https://www.biostars.org/p/129763/)
```
samtools bam2fq <BAM> | seqtk seq -A > <FASTA>
```

###[BAM to BED] (https://www.biostars.org/p/85990/)
```
bedtools bamtobed -i data.bam > data.bed 
samtools view data.bam > data.sam
paste data.bed data.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' | head
```
###[Check if BED file is tab-delimited] (https://www.biostars.org/p/127275/)
```
sed 's/ \+/\t/g'
```

<h1><p align="center">Analyze Alignment Data</p></h1>
###[Calculate coverage with "genomecov"] (http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
* [Calculate coverage of 5' end of reads aligned to genome] (https://www.biostars.org/p/80236/)
```
bedtools genomecov [-ibam] <BAM> -g <GENOME>
```

###[Obtain upstream and downstream sequences from a FASTA file] (https://www.biostars.org/p/82776/)
* [BEDtools flank] (http://bedtools.readthedocs.io/en/latest/content/tools/flank.html)
* [BEDtools getfasta] (http://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)

<h1><p align="center">Additional Helpful Resources</p></h1>

* [Reference for I/O redirection] (http://www.tldp.org/LDP/abs/html/io-redirection.html)
* [Extract filename from filepath] (http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash)
* [Extract directory from filepath] (http://stackoverflow.com/questions/6121091/get-file-directory-path-from-filepath)
* [Sequencing Contamination] (http://seqanswers.com/forums/showthread.php?t=12520)
* [Duplicate reads vs. Singletons] (http://sfg.stanford.edu/quality.html)
* Zero-based vs. One-based files
*  https://www.biostars.org/p/51504/
*  https://www.biostars.org/p/84686/

