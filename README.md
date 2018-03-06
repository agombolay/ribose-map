# Ribose-Map Bioinformatics Toolkit
## Toolkit for profiling the identity and distribution of rNMPs embedded in DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

## Modules:
1. **Alignment**: Align reads to the reference with Bowtie2 and de-depulicated based on UMI's UMI-tools
2. **Coordinates**: Locate genomic coordinates of rNMPs for ribose-seq, Pu-seq, emRibo-seq, or HydEn-seq
3. **Frequencies**: Calculate and visualize frequencies of nucleotides at and flanking sites of rNMP incorporation
4. **Distribution**: Visualize coverage of rNMPs across chromosomes and create bedgraph files for genome browser

## How to set up repository:
git clone https://github.com/agombolay/ribose-map/

* It is recommended to add the scripts to your $PATH  
* Mitochondria should be named chrM or MT in FASTA 

## Required software dependencies:
* [Bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.1), [UMI-tools](https://github.com/CGATOxford/UMI-tools), [SAMtools](http://www.htslib.org/download/), [BEDtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), [cutadapt](http://cutadapt.readthedocs.io/en/stable/), [Python](https://www.python.org/), and [R](https://cran.r-project.org/) (required libraries: tools, optparse and ggplot2)

## Command usage:

|   Alignment Module   |   Coordinates Module   |   Frequencies Module   |   Distribution Module   |
| -------------------- | ---------------------- | ---------------------  | ----------------------- |
| Alignment.sh config  | Coordinates.sh config  | Frequencies.sh config  | Distribution.sh config  |
|                      |                        | Frequencies.R config   | Distribution.R config   |

## Example config file:
```
#Optional
barcode='TCA'
umi='NNNNNNXXXNN'

#Required
sample='sample1'
technique='ribose-seq'
directory='/filepath/ribose-map'
index='/filepath/ribose-map/indices/sacCer2'
reference='/filepath/ribose-map/references/sacCer2.fa'
read1='/filepath/ribose-map/trimmed/FS120-trimmed1.fastq'
read2='/filepath/ribose-map/trimmed/FS120-trimmed2.fastq'
```
