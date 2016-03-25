# Ribose-seq-Project
Alli Gombolay, M.P.H  
Storici Lab | School of Biology  
Georgia Institute of Technology  

##Project Overview
###Non-LSF Dependent Version of Ribose-seq Analysis Pipeline  

**References**:  
[Ribose-seq *Nature Methods* Paper, 2015]
(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4686381/pdf/nihms742750.pdf)  
[Georgia Tech News Article on Ribose-seq]
(http://www.news.gatech.edu/2015/01/26/ribose-seq-identifies-and-locates-ribonucleotides-genomic-dna)

##Required Input Files:  
-Sequencing data in FASTQ file format  
-Yeast data: http://amc-sandbox.ucdenver.edu/User13/outbox/2016/  

##Software Requirements:  
* [Bowtie] (https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/): Align raw sequencing reads to reference genome
 * [How to test if pre-built Bowtie index is properly installed] (http://bowtie-bio.sourceforge.net/tutorial.shtml)
* [bedToBigBed and bedGraphToBigWig] (http://hgdownload.cse.ucsc.edu/admin/exe/)
* [bedtools]  (http://bedtools.readthedocs.org/en/latest/content/installation.html)
* [umitools] (https://github.com/brwnj/umitools)
* [SAMtools] (http://www.htslib.org/download/)
* [Python] (https://www.python.org/downloads/)  
* [R]  (https://www.r-project.org/)

##Set-up:
###Part A: Software Set-up  
-Script to download and install software automatically  

###Part B: Directory Set-up  
1. Clone the Ribose-seq Analysis Pipeline Directory Structure:  
```git clone https://github.com/agombolay/Ribose-seq-Project/tree/master/Ribose-seq-Directory```
