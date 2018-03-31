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

1. **Download Git repository**:  
    Click [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) for information on installing Git
    ```
    git clone https://github.com/agombolay/ribose-map/
    ```
  * It is recommended to add the scripts to your $PATH  
  * Mitochondria should be named chrM or MT in FASTA 

  Ribose-Map uses [Bowtie 2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.1), [BEDtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), [SAMtools](http://www.htslib.org/download/), [cutadapt](http://cutadapt.readthedocs.io/en/stable/), [UMI-tools](https://github.com/CGATOxford/UMI-tools), [R](https://cran.r-project.org/), and [Python](https://www.python.org/) to analyze and visualize data.  
To ensure easy installation and versioning of this software, we recommend using the MiniConda package manager.

2. **Install pre-requisites for conda**:
     ```
     python3 -m pip install pycosat pyyaml requests --user
     ```

3. **Install MiniConda and download third party software**:

     1. Install MiniConda:  
        Download .sh script from https://conda.io/miniconda.html  
        Follow prompts to install Miniconda and add it to your PATH 
        ```
        sh Miniconda3-latest-Linux-x86_64.sh
        ```
        
     2. Source your .bashrc to ensure that MiniConda loads properly:
        ```
        source ~/.bashrc
        ```

     3. Verify Miniconda was installed (Installed packages are displayed):
        ```
        conda list
        ```

     4. Update conda and install anaconda-client:
        ```
        conda update conda
        conda install anaconda-client anaconda-build conda-build
        ```
        
      5. Create conda environment for Ribose-Map:  
         ribosemap_env.yaml is located in ribose-map/lib
         ```
         conda env create -n ribosemap_env --file ribosemap_env.yaml
         ```

4. **Activate conda environment to use Ribose-Map**:
```
source activate ribosemap_env
```

5. **Once the analysis is complete, exit environment**:  
```
source deactivate ribosemap_env
```

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
