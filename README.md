<img src="https://github.com/agombolay/ribose-map/blob/master/logo.pdf" width="400px" height="300px" />

# Ribose-Map Bioinformatics Toolkit
## Toolkit for mapping rNMPs embedded in DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

## Modules:
**Processing rNMP sequencing data**:
* **Alignment**: Align reads to the reference with Bowtie2 and de-depulicated based on UMI's UMI-tools
* **Coordinate**: Locate genomic coordinates of rNMPs for ribose-seq, Pu-seq, emRibo-seq, or HydEn-seq

**Analyzing genome-wide distribution of rNMPs and their sequence context**:
* **Sequence**: Calculate and visualize frequencies of nucleotides at and flanking sites of rNMP incorporation
* **Distribution**: Visualize coverage of rNMPs across genome and create bedgraph files for genome browser

## How to set up repository:
```
git clone https://github.com/agombolay/ribose-map/
```

* It is recommended to add the scripts to your $PATH  
* Mitochondria should be named chrM or MT in FASTA 

## Software dependencies:
### Required software:
* [Bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.1), [BEDtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), [SAMtools](http://www.htslib.org/download/), and [R](https://cran.r-project.org/) (tools and ggplot2)

### Additional software:
* [cutadapt](http://cutadapt.readthedocs.io/en/stable/) is required if libraries contain a 5' molecular barcode
* [UMI-tools](https://github.com/CGATOxford/UMI-tools) is required if libraries contain a unique molecular identifier
  * Note: UMI-tools and cutadapt both require [Python](https://www.python.org/) to install and run

## Command usage:

| Alignment Module        | Coordinate Module       | Sequence Module         | Distribution Module     |
| ----------------------- | ----------------------- | ----------------------- | ----------------------- |
| alignment.sh config     | coordinate.sh config    | sequence.sh config      | distribution.sh onfig   |
|                         |                         | sequence.R config       | distribution.R config   |

## Example config:
```
#Sample name
sample='sample1'

#Library prep
barcode='TCA'
pattern='NNNNNNXXXNN'

#rNMP Sequencing
technique='ribose-seq'

#Reference genome
fasta='/filepath/sacCer2.fa'
basename='/filepath/sacCer2'

#FASTQs files of reads
read1='/filepath/sample1_1.fastq'
read2='/filepath/sample1_2.fastq'

#Ribose-Map repository
repository='/filepath/ribose-map'
```
