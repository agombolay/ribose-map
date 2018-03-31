![Logo](https://github.com/agombolay/ribose-map/blob/master/logo.png)
# A bioinformatics toolkit for mapping rNMPs embedded in DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

## Modules:
**Processing rNMP sequencing data**:
* **Alignment**: Align reads to the reference with Bowtie2 and de-depulicated based on UMI's UMI-tools
* **Coordinate**: Locate genomic coordinates of rNMPs for ribose-seq, Pu-seq, emRibo-seq, or HydEn-seq

**Analyzing genome-wide distribution of rNMPs and their sequence context**:
* **Sequence**: Calculate and visualize frequencies of nucleotides at and flanking sites of embedded rNMPs
* **Distribution**: Visualize coverage of rNMPs across genome and create bedgraph files for genome browser

## How to set up repository:

Download git repository:
```
git clone https://github.com/agombolay/ribose-map/
```
* It is recommended to add the scripts to your $PATH  
* Mitochondria should be named chrM or MT in FASTA 

Ribose-Map uses several standard bioinformatics tools for data analysis and R's ggplot2 for visualizing the results.  
To ensure easy installation and versioning of this software, we recommend using the MiniConda package manager.

Install conda prerequisites:
```
python3 -m pip install pycosat pyyaml requests --user
```

Install MiniConda:
```
sh Miniconda3-latest-Linux-x86_64.sh
```

Press ENTER when prompted; when asked for installation path, type yes and press ENTER to use your HOME folder as the site of installation of Miniconda3 or enter path to the folder where you want it to be installed.

When asked if you want to add Miniconda3 to your .bashrc, type yes and press ENTER.

Source your .bashrc:
```
source ~/.bashrc
```

Verify Miniconda was installed:
Installed packages are displayed.
```
conda list
```

Update conda:
```
conda update conda
```

Install anaconda client:
```
conda install anaconda-client anaconda-build conda-build
```

Create conda environment for Ribose-Map:
```
conda env create -n ribosemap_env --file ribosemap_env.yaml
```

Activate conda environment to use Ribose-Map:
```
source activate ribosemap_env
```

Once the analysis is complete, exit the environment:  
```
source deactivate ribosemap_env
```

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
