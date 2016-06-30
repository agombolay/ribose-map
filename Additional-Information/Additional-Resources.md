##Additional Helpful Resources:  

* [Reference for I/O redirection] (http://www.tldp.org/LDP/abs/html/io-redirection.html)
* [Extract filename from filepath] (http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash)
* [Extract directory from filepath] (http://stackoverflow.com/questions/6121091/get-file-directory-path-from-filepath)

###[Duplicate reads vs. Singletons] (http://sfg.stanford.edu/quality.html)

###Zero-based vs. One-based files
* https://www.biostars.org/p/51504/
* https://www.biostars.org/p/84686/

###[Make sure BED file is tab-delimited] (https://www.biostars.org/p/127275/)
```
sed 's/ \+/\t/g'
```

###[bam-readcount] (https://github.com/genome/bam-readcount)

###[Obtain upstream and downstream sequences from a FASTA file] (https://www.biostars.org/p/82776/)
* [BEDtools flank] (http://bedtools.readthedocs.io/en/latest/content/tools/flank.html)
* [BEDtools getfasta] (http://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html)

####[Calculate coverage with "genomecov"] (http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
* [Calculate coverage of 5' end of reads aligned to genome] (https://www.biostars.org/p/80236/)
```
bedtools genomecov [-ibam] <BAM> -g <GENOME>
```

##[Sequencing Contamination] (http://seqanswers.com/forums/showthread.php?t=12520)

####[BAM to FASTA:] (https://www.biostars.org/p/129763/)
```
samtools bam2fq <BAM> | seqtk seq -A > <FASTA>
```

* [Count number of reads in a FASTA file] (http://thegenomefactory.blogspot.com/2011/09/counting-sequences-with-unix-tools.html)

##Check quality of reads
####[FastQC] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

##Trim reads based on quality and adapters
Importance of removing adapters: Increase mapping percentage and decrease incorrect mappings

####[Trimmomatic] (http://www.usadellab.org/cms/?page=trimmomatic)
* [Manual for Trimmomatic] (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

##[Count number of reads in BAM file] (http://crazyhottommy.blogspot.com/2013/06/count-how-many-mapped-reads-in-bam-file.html)

**Mapped reads**:
```
samtools view -c -F 4
```

**Unmapped reads**:
```
samtools view -c -f 4
```

##View alignment data

####[SAMtools tview] (http://samtools.sourceforge.net/tview.shtml)

```
samtools tview <BAM> --reference <Fasta>
```
