#Analysis of Ribose-seq Libraries
Alli Gombolay, M.P.H  
[Storici Lab] (http://www.storicilab.gatech.edu/) | [School of Biology] (http://www.biology.gatech.edu/)  
[Georgia Institute of Technology] (http://www.gatech.edu/)  
Contact: alli.gombolay@gatech.edu

Last Updated: May 2016  

##Project Overview
**Goal**: To create a streamlined bioinformatics pipeline to analyze Ribose-seq libraries on any Linux platform

**Background**: The original scripts used to analyze Ribose-seq libraries required LSF job scheduling software to run properly.  To enable the user to efficiently analyze their input sequenced Ribose-seq libraries, I have created a streamlined, non-LSF dependent version of the Ribose-seq analysis pipeline.  This pipeline will run on any Linux platform.  I have also included all required input files, the appropriate directory structure that can be cloned as a GitHub repository, and a script to automatically download and install the required software.  I have tailored the pipeline for four different species, including *Saccharomyces cerevisiae*, *Escherichia coli K12*, *Mus musculus*, and humans.

**Significance**: This pipeline allows any user to easily download and install all required software, download all required input files, and download the appropriate directory structure to analyze Ribose-seq libraries.  Now, users can readily choose the species of interest and run the corresponding pipeline to map the distribution of rNMPs in their Ribose-seq libraries.  

**Applications**:  

##**References**  
* [Ribose-seq *Nature Methods* Paper, 2015]
(http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4686381/pdf/nihms742750.pdf)  
* [Jay Hesselberth's GitHub Page for Ribose-seq]
(https://github.com/hesselberthlab/modmap/tree/snake/pipeline/ribose-seq-ms)
* [Georgia Tech 2015 News Article on Ribose-seq]
(http://www.news.gatech.edu/2015/01/26/ribose-seq-identifies-and-locates-ribonucleotides-genomic-dna)

##Program Set-up
###Part A: Software Set-up  
* Script to download and install software automatically
* [Required Software] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Required-Software.md)
* [Required Input Files] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Required-Input-Files.md)

###Part B: Directory Set-up  
1. Clone the Ribose-seq Analysis Pipeline Directory Structure:  
```git clone https://github.com/agombolay/Ribose-seq-Project/tree/master/Ribose-seq-Directory```
