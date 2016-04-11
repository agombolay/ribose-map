# Ribose-seq-Project
Alli Gombolay, M.P.H  
Storici Lab | School of Biology  
Georgia Institute of Technology  

##Project Overview
###Non-LSF Dependent Version of Ribose-seq Analysis Pipeline  

**References**:  
* [Ribose-seq *Nature Methods* Paper, 2015]
(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4686381/pdf/nihms742750.pdf)  
* [Georgia Tech News Article on Ribose-seq]
(http://www.news.gatech.edu/2015/01/26/ribose-seq-identifies-and-locates-ribonucleotides-genomic-dna)
* [Jay Hesselberth's GitHub Page]
(https://github.com/hesselberthlab/modmap/tree/snake/pipeline/ribose-seq-ms)

##Biological Information on *Saccharomyces cerevisiae*
* Chromosomes: 16
* Genome size: 12 million base pairs  
(Reference: https://www.ncbi.nlm.nih.gov/genome/15)

##Required Input Files:  
-Sequencing data in FASTQ file format  
-Saccharomyces cerevisiae data: http://amc-sandbox.ucdenver.edu/User13/outbox/2016/  

##Convert input raw sequencing files from SRA to FASTQ format
* SRA Toolkit (fastq-dump):  
http://www.ncbi.nlm.nih.gov/books/NBK158900/  
http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software  
http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc  

##Convert input reference genome files from 2bit to fasta format
* http://hgdownload.soe.ucsc.edu/admin/exe/

##Set-up:
###Part A: Software Set-up  
-Script to download and install software automatically  

###Part B: Directory Set-up  
1. Clone the Ribose-seq Analysis Pipeline Directory Structure:  
```git clone https://github.com/agombolay/Ribose-seq-Project/tree/master/Ribose-seq-Directory```
