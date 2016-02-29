# Ribose-seq-Project
Alli Gombolay, M.P.H  
Storici Lab | School of Biology  
Georgia Institute of Technology  

##Project Overview
###Open-source Version of Ribose-seq Analysis Pipeline  

**References**:  
[Ribose-seq *Nature Methods* Paper, 2015]
(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4686381/pdf/nihms742750.pdf)  
[Georgia Tech News Article on Ribose-seq]
(http://www.news.gatech.edu/2015/01/26/ribose-seq-identifies-and-locates-ribonucleotides-genomic-dna)

##Required Files:  
-Input Files: FASTQ files  
-http://amc-sandbox.ucdenver.edu/User13/outbox/2016/  

##Software Requirements:  
-bowtie (https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/)  
-umitools (https://github.com/brwnj/umitools); Also requires pysam and editdist  
-bedtools  (http://bedtools.readthedocs.org/en/latest/content/installation.html)  
-bedToBigBed and bedGraphToBigWig (http://hgdownload.cse.ucsc.edu/admin/exe/)  
-samtools (http://www.htslib.org/download/)  
-python (https://www.python.org/downloads/)  
-R  (https://www.r-project.org/)  

##Getting ready:
Part A: Software Set-up
1. Create a bin folder in home directory:  
```mkdir bin```  

2. If not already done, add bin folder to PATH:  
```echo "export PATH="~/bin:$PATH"" >> ~/.bashrc```  

3. Download all required software into bin folder:  
```cd path/to/bin```  

Part B: Directory Set-up:  
1. Clone the Ribose-seq Analysis Pipeline Directory Structure:  
```git clone https://github.com/agombolay/Ribose-seq-Project/tree/master/Ribose-seq-Directory```  

##Order of Execution of Scripts:  
1.  
2.  
3.  
4.  
