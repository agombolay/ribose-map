<h2><p align="center">Count number of reads</p></h2>
### [Count number of reads in FASTA file](http://thegenomefactory.blogspot.com/2011/09/counting-sequences-with-unix-tools.html)

```
grep '>' in.fasta | wc -l
```

### [Count number of reads in BAM file](http://crazyhottommy.blogspot.com/2013/06/count-how-many-mapped-reads-in-bam-file.html)

**Mapped reads**:
```
samtools view -c -F 4
```

**Unmapped reads**:
```
samtools view -c -f 4
```

<h2><p align="center">Check quality of reads</p></h2>
### [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

* [Example Reports (Good, bad, adapter contaminated, etc.)]
(http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Information reagarding overrepresented Sequences] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html)

<h2><p align="center">Trim reads based on quality and adapters</p></h2>
**Importance of removing adapters**: Increase mapping percentage and decrease incorrect mappings

### [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  * [Manual for Trimmomatic] (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

```
java -jar trimmomatic-0.36.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
* ILLUMINACLIP:TruSeq3-SE.fa:2:30:10: Trims adapters from reads
* TRAILING:3: Cuts bases off the end of read, if below threshold quality
 * Do not use LEADING:3 opton since the UMI need to be retained
* MINLEN:36: Drop the read if it is below a specified length

<h2><p align="center">View alignment data</p></h2>

### [SAMtools tview](http://samtools.sourceforge.net/tview.shtml)

```
samtools tview -p <chr:pos> <BAM> --reference <Fasta>
```
* ["."](http://samtools.sourceforge.net/pileup.shtml) = match to the reference base on the forward strand
* [","](http://samtools.sourceforge.net/pileup.shtml) = match to the reference base on the reverse strand
* "?" = help menu

<h2><p align="center">Convert file formats</p></h2>

### [BAM to FASTA](https://www.biostars.org/p/129763/)
```
samtools bam2fq <BAM> > <output FASTQ>
seqtk seq -A <input FASTQ> > <output FASTA>
```

### [BAM to BED](https://www.biostars.org/p/85990/)
```
bedtools bamtobed -i data.bam > data.bed 
samtools view data.bam > data.sam
paste data.bed data.sam | awk -v "OFS=\t" '{print $1, $2, $3, $16, $6}' | head
```
### [Check if BED file is tab-delimited](https://www.biostars.org/p/127275/)
```
sed 's/ \+/\t/g'
```

<h2><p align="center">Analyze alignment data</p></h2>
### [Calculate coverage with "genomecov"](http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
* [Calculate coverage of 5' end of reads aligned to genome] (https://www.biostars.org/p/80236/)
```
bedtools genomecov -d 5 [-ibam] <BAM> -g <GENOME>
```

### [Obtain upstream and downstream sequences from a FASTA file](https://www.biostars.org/p/82776/)
* [BEDtools flank] (http://bedtools.readthedocs.io/en/latest/content/tools/flank.html)
```
bedtools flank -i <BED> -g <GENOME> -b 100
```

* [BEDtools getfasta](http://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)
```
bedtools getfasta -fi <input FASTA> -bed <BED> -fo <output FASTA>
```
* Chromosome name in header of FASTA file must match those in BED file

<h2><p align="center">Additional helpful resources</p></h2>

### Programming
* [Reference for I/O redirection] (http://www.tldp.org/LDP/abs/html/io-redirection.html)
* [Extract filename from filepath] (http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash)
* [Extract directory from filepath] (http://stackoverflow.com/questions/6121091/get-file-directory-path-from-filepath)

### Bioinformatics
* [Sequencing Contamination] (http://seqanswers.com/forums/showthread.php?t=12520)
* [Duplicate reads vs. Singletons] (http://sfg.stanford.edu/quality.html)
* Zero-based vs. One-based files:
 * https://www.biostars.org/p/51504/
 * https://www.biostars.org/p/84686/
