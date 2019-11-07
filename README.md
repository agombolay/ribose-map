![Logo](https://github.com/agombolay/Images/blob/master/Logo.png)
# A bioinformatics toolkit for mapping rNMPs in genomic DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

**If you use Ribose-Map, please use the following citation**:  
Gombolay, AL, FO Vannberg, and F Storici. Ribose-Map: a bioinformatics toolkit to map ribonucleotides embedded in genomic DNA. *Nucleic Acids Research* 2019. DOI: 10.1093/nar/gky874.

Ribose-Map is the first and only bioinformatics toolkit that can process and analyze high-throughput rNMP sequencing data that was generated using any of the currently available rNMP sequencing techniques. Current techniques include: Alk-HydEn-seq, emRiboSeq, Pu-seq, RADAR-seq, ribose-seq, RHII-HydEn-seq.

## Modules
**Alignment**: Aligns reads to reference genome, de-depulicates based on UMI, and de-multiplexes by barcode  
**Coordinate**: Calculates single-nucleotide chromosomal coordinates of rNMPs based on sequencing technique used  

Based on the chromosomal coordinates, the Composition, Sequence, and Distribution Modules...  
**Composition**: Calculates and plots the percentage of r{A,C,G,U}MP noramlized to the reference genome sequence  
**Sequence**: Calculates and plots the nucleotide frequencies of rNMP sites and those up/downstream from those sites  
**Distribution**: Creates bedgraph files that can be uploaded to a genome browser and plots per-nucleotide rNMP coverage
 
&nbsp;
## Required Data
**Ribose-Map requires the following files**:
1. FASTQ files of sequencing data generated using NGS
2. Bowtie2 indexes (extension .bt2) for reference genome
3. FASTA file of nucleotide sequence of reference genome
4. Configuration file (example is provided in ribose-map/lib)

Note: If you have a list of single nucleotide coordinates that were not created using the Coordinate Module and would like to use the Sequence and/or Distribution Modules, you will need to first create a directory with the path shown below:

"$repository/results/$sample/coordinate" ($repository and $sample variables should be those provided in the config file)

&nbsp;
## Software Installation

1. **Clone Ribose-Map GitHub repository**:  
   * Click [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) for information on installing Git on Linux
   ```
   git clone https://github.com/agombolay/ribose-map/
   ```
   
    Ribose-Map uses [Bowtie 2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.1), [BEDtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), [SAMtools](http://www.htslib.org/download/), [cutadapt](http://cutadapt.readthedocs.io/en/stable/), [UMI-tools](https://github.com/CGATOxford/UMI-tools), [R](https://cran.r-project.org/), and [Python](https://www.python.org/) to analyze and visualize data.  
To ensure easy installation and versioning of these software, we recommend using the MiniConda package manager.

2. **Install pre-requisite packages for conda**:
   ```
   python3 -m pip install pycosat pyyaml requests --user
   ```

3. **Install MiniConda and source your .bashrc file**:  
   * Follow prompts to install MiniConda in home directory
   * Run these commands from within ribose-map/lib directory
   ```
   sh Miniconda3-latest-Linux-x86_64.sh && source ~/.bashrc
   ```

4. **Create environment in which to run Ribose-Map**:  
   * Run these commands from within ribose-map/lib directory
   ```
   conda update conda
   conda install anaconda-client anaconda-build conda-build
   conda env create --name ribosemap_env --file ribosemap.yaml
   ```

&nbsp;
## How to run Ribose-Map
Before proceeding, close current terminal window and open a new window to refresh the settings  
* If SAMtools gives an error, then install this dependency: conda install -c conda-forge ncurses

1. **Activate environment to access software**:
```
source activate ribosemap_env
```

2. **Run scripts with configuration_file as input**:
```
ribosemap {alignment, coordinate, composition, sequence, distribution} config
```

3. **Once the analysis is complete, exit environment**:  
```
source deactivate ribosemap_env
```
