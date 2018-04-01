![Logo](https://github.com/agombolay/ribose-map/blob/master/logo.png)
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

1. **Download Git repository**:  
   * Click [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) for information on installing Git on Linux
   * Add the /ribose-map/modules directory to your PATH
    ```
    git clone https://github.com/agombolay/ribose-map/
    ```

    Ribose-Map uses [Bowtie 2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.1), [BEDtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), [SAMtools](http://www.htslib.org/download/), [cutadapt](http://cutadapt.readthedocs.io/en/stable/), [UMI-tools](https://github.com/CGATOxford/UMI-tools), [R](https://cran.r-project.org/), and [Python](https://www.python.org/) to analyze and visualize data.  
To ensure easy installation and versioning of this software, we recommend using the MiniConda package manager.

2. **Install pre-requisites for conda**:
     ```
     python3 -m pip install pycosat pyyaml requests --user
     ```

3. **Install MiniConda and software dependencies**:  
   * Example configuration_file is available at lib/config
     1. Install MiniConda and source .bashrc:  
        * Follow series of prompts (press ENTER and type yes)
        ```
        sh lib/Miniconda3-latest-Linux-x86_64.sh && source ~/.bashrc
        ```

     2. Create conda environment for Ribose-Map:  
        * Software dependencies will be installed in environment
        ```
        conda update conda
        conda install anaconda-client anaconda-build conda-build
        conda env create -n lib/ribosemap_env --file ribosemap_env.yaml
        ```

&nbsp;
## How to run Ribose-Map
1. **Activate conda environment to use Ribose-Map**:
```
source activate ribosemap_env
```

2. **Run .sh and .R scripts along with configuration_file**:
   * Example configuration_file is available at lib/config

     1. **Process the data**
     ```
     alignment.sh configuration_file
     coordinate.sh configuration_file
     ```
     2. **Analyze the data**
     ```
     sequence.sh configuration_file
     distribution.sh configuration_file
     ```
     3. **Visualize the results**
     ```
     sequence.R configuration_file
     distribution.R configuration_file
     ```

3. **Once the analysis is complete, exit environment**:  
```
source deactivate ribosemap_env
```
