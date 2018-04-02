![Logo](https://github.com/agombolay/Images/blob/master/logo.png)
# A bioinformatics toolkit for mapping rNMPs embedded in DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

## Modules
**Process rNMP sequencing data independent of the rNMP sequencing technique**:  
* **Alignment**: Align reads to the reference with Bowtie2 and de-depulicated based on UMI's UMI-tools  
* **Coordinate**: Locate genomic coordinates of rNMPs for ribose-seq, Pu-seq, emRibo-seq, or HydEn-seq  

**Analyze sequence context of embedded rNMPs and their genome-wide distribution**:  
* **Sequence**: Calculate and visualize frequencies of nucleotides at and flanking sites of embedded rNMPs  
* **Distribution**: Visualize coverage of rNMPs across genome and create bedgraph files for genome browser  
 
&nbsp;
## Required Data
**Ribose-Map requires the following files**:
1. FASTQ files of sequencing data generating using NGS
2. Bowtie2 indexes (extension .bt2) for reference genome
3. FASTA file of nucleotide sequence of reference genome

&nbsp;
## Software Installation

1. **Clone Ribose-Map GitHub repository**:  
   * Click [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) for information on installing Git on Linux
   * Add the ribose-map/modules directory to your PATH
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
   ```
   sh lib/Miniconda3-latest-Linux-x86_64.sh && . ~/.bashrc
   ```

4. **Create environment in which to run Ribose-Map**:  
   ```
   conda update conda
   conda install anaconda-client anaconda-build conda-build
   conda env create -n ribosemap_env --file lib/ribosemap.yaml
   ```

&nbsp;
## How to run Ribose-Map
1. **Activate environment to access software**:
```
source activate ribosemap_env
```

2. **Run .sh and .R scripts along with configuration_file**:
   * An example configuration file is available at lib/configuration_file
   * .sh and R scripts should be run within ribose-map/modules directory

     1. **Process the data**
     ```
     sh alignment.sh configuration_file
     sh coordinate.sh configuration_file
     ```
     2. **Analyze the data**
     ```
     sh sequence.sh configuration_file
     sh distribution.sh configuration_file
     ```
     3. **Visualize the results**
     ```
     R sequence.R configuration_file
     R distribution.R configuration_file
     ```

3. **Once the analysis is complete, exit environment**:  
```
source deactivate ribosemap_env
```
