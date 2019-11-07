![Logo](https://github.com/agombolay/Images/blob/master/Logo.png)
# A bioinformatics toolkit for mapping rNMPs in genomic DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

**If you use Ribose-Map, please use the following citation**:  
Gombolay, AL, FO Vannberg, and F Storici. Ribose-Map: a bioinformatics toolkit to map ribonucleotides embedded in genomic DNA. *Nucleic Acids Research*, Volume 47, Issue 1, 10 Jan 2019, Page e5, https://doi.org/10.1093/nar/gky874

Ribose-Map is the first and only bioinformatics toolkit that can process and analyze high-throughput rNMP sequencing data that was generated using any of the currently available rNMP sequencing techniques. Current techniques include: Alk-HydEn-seq, emRiboSeq, Pu-seq, RADAR-seq, ribose-seq, RHII-HydEn-seq.

In addition to analyzing rNMP sequencing data, Ribose-Map's Composition, Sequence, Distribution, and Hotspot Modules can also be used to analyze any type of single-nucleotide sequencing data.  

## Modules
The **Alignment Module** aligns reads to the reference genome and de-depulicates/de-multiplexes the reads

The **Coordinate Module** calculates single-nucleotide chromosomal coordinates of rNMPs based on the positions of the aligned reads. The coordinates can be calculated for any of the currently available rNMP sequencing techniques.  
  
Based on the chromosomal coordinates of rNMPs, the **Composition**, **Sequence**, **Distribution**, and **Hotspot Modules**...  

* Calculate and plot the percentage of r{A,C,G,U}MP normalized to the reference genome sequence  

* Calculate and plot the nucleotide frequencies of rNMP sites and up/down-stream from those sites  

* Create bedgraph files of data to be uploaded and visualized in any genome browser of choice and also plot the per-nucleotide rNMP coverage for each chromosome  

* Calculate the top 1% of rNMP sites and creates sequence logo plots for those sites

## Required Data
1. FASTQ file of rNMP sequencing data
2. FASTA file of nucleotide sequence of reference genome
3. Bowtie2 indexes (extension = .bt2) of reference genome
4. Configuration file (example provided in ribose-map/lib)

**Note**: If you have a list of single nucleotide coordinates that were not created using the Coordinate Module and would like to directly use the Sequence, Distribution, or Hotspot Modules, first create a directory with the path shown below:

"$repository/results/$sample/coordinate$quality" ($variables should be those provided in config)

## Install Software

1. **Clone Ribose-Map repository**:  
   ```
   git clone https://github.com/agombolay/ribose-map/
   ```
   
    Ribose-Map uses [Bowtie 2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.1), [BEDtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), [SAMtools](http://www.htslib.org/download/), [cutadapt](http://cutadapt.readthedocs.io/en/stable/), [UMI-tools](https://github.com/CGATOxford/UMI-tools), [R](https://cran.r-project.org/), and [Python](https://www.python.org/) to analyze and visualize data.  
To ensure easy installation and versioning of these software, we recommend using MiniConda package manager.

2. **Install pre-requisite packages for conda**:
   ```
   python3 -m pip install pycosat pyyaml requests --user
   ```

3. **Install MiniConda and source .bashrc file**:  
   ```
   sh Miniconda3-latest-Linux-x86_64.sh && source ~/.bashrc
   ```

4. **Create conda environment**:  
   ```
   conda update conda
   conda install anaconda-client anaconda-build conda-build
   conda env create --name ribosemap_env --file ribosemap.yaml
   ```

## Run Ribose-Map
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
