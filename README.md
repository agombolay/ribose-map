# Computational Analysis of Ribose-seq Libraries
Alli Gombolay, M.P.H  
[Storici Lab](http://www.storicilab.gatech.edu/) | [School of Biology](http://www.biology.gatech.edu/)  
[Georgia Institute of Technology](http://www.gatech.edu/)  
Contact: alli.gombolay@gatech.edu

Last Updated: May 2016  

## Project Overview
**Goal**: To create a user-friendly bioinformatics toolkit to analyze ribose-seq libraries

<p align="justify">
<b>Background</b>: The original scripts used to analyze Ribose-seq libraries required Platform Load Sharing Facility (LSF) workload management software to run properly.  To enable the user to efficiently analyze their sequenced ribose-seq libraries, I have created an updated, streamlined version of the Ribose-seq analysis pipeline that does not depend on LSF.  This pipeline has been designed to run on any Linux platform, and it consists of both Bash scripts and R scripts.  The pipeline includes a script to automatically download all required input files, the appropriate directory structure (as a GitHub repository), and the required software.
</p>

<p align="justify">
<b>Significance</b>: With minimal software dependencies, environment set-up, and comprehensive documentation, Ribose-Map allows users to readily profile the identity and distribution of rNMPs in any organism of interest.
</p>

## Usage
To display help menu: `script -h`  

## References  
* [Ribose-seq *Nature Methods* Paper, 2015](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4686381/pdf/nihms742750.pdf)  
* [Jay Hesselberth's GitHub Page for Ribose-seq](https://github.com/hesselberthlab/modmap/tree/snake/pipeline/ribose-seq-ms)
* [Georgia Tech 2015 News Article on Ribose-seq](http://www.news.gatech.edu/2015/01/26/ribose-seq-identifies-and-locates-ribonucleotides-genomic-dna)

## Program Set-up
### Part A: Software Set-up  
* Script to automatically download required input files
* Script to automatically download and install software

### Part B: Directory Set-up  
1. Clone the Ribose-seq Analysis Pipeline Directory Structure:  
```git clone https://github.com/agombolay/Ribose-seq-Project/tree/master/Ribose-seq-Directory```
