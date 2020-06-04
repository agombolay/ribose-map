![Logo](https://github.com/agombolay/Images/blob/master/Logo.png)
# A bioinformatics toolkit for mapping rNMPs in genomic DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

**If you use Ribose-Map, please use the following citation**:  
Gombolay, AL, FO Vannberg, and F Storici. Ribose-Map: a bioinformatics toolkit to map ribonucleotides embedded in genomic DNA. *Nucleic Acids Research*, Volume 47, Issue 1, 10 Jan 2019, Page e5, https://doi.org/10.1093/nar/gky874. 

## Modules
The **Alignment Module** aligns reads to reference genome and de-depulicates/de-multiplexes reads if needed
usage: ribose-map/modules/ribosemap alignment config

The **Coordinate Module** calculates single-nucleotide chromosomal coordinates of rNMPs based on aligned reads for any NGS rNMP sequencing technique, including ribose-seq, emRiboSeq, RHII-HydEn-seq, Alk-HydEn-seq, or Pu-seq.  
usage: ribose-map/modules/ribosemap coordinate config

Based on the coordinates of rNMPs, ...  

* **Composition**: Plots the percentage of r{A,C,G,U} normalized to the reference genome sequence  
usage: ribose-map/modules/ribosemap composition config

* **Sequence**: Plots the nucleotide frequencies of rNMP sites and up/down-stream from those sites  
usage: ribose-map/modules/ribosemap sequence config

* **Distribution**: Creates bedgraph files of per-nucleotide rNMP coverage to be visualized in a genome browser and plots the per-nucleotide rNMP coverage for each chromosome (%) normalized to account for sequencing depth
usage: ribose-map/modules/ribosemap distribution config

* **Hotspot**: Calculates the top 1% most abundant sites of rNMP incorporation and creates sequence logos
usage: ribose-map/modules/ribosemap hotspot config

## Required Data
1. FASTQ file of rNMP sequencing data
2. FASTA file of nucleotide sequence of reference genome
3. Bowtie2 indexes (extension = .bt2) of reference genome
4. Configuration file (example provided in ribose-map/lib)

**Note**: If you have a list of single nucleotide coordinates that were not created using the Coordinate Module and would like to directly use the Sequence, Distribution, or Hotspot Modules, first create a directory with the path shown below:

"$repository/results/$sample/coordinate$quality" ($variables should be those provided in config)

## Install Software

1. **Create conda software environment**:  
* If you do not already have MiniConda installed, please visit [this link](https://docs.conda.io/en/latest/miniconda.html)
   ```bash
   conda update conda
   conda install anaconda-client anaconda-build conda-build
   conda env create --name ribosemap_env --file ribosemap.yaml
   ```

2. **Clone Ribose-Map GitHub repository**:  
   ```bash
   git clone https://github.com/agombolay/ribose-map/
   ```
   
## Run Ribose-Map
Before proceeding, close current terminal window and open a new window to refresh the settings  
* If SAMtools gives an error, then install this dependency: conda install -c conda-forge ncurses

1. **Activate environment to access software**:
```bash
source activate ribosemap_env
```

2. **Run scripts with configuration_file as input**:
```bash
ribosemap {alignment, coordinate, composition, sequence, distribution, hotspot} config
```

3. **Once the analysis is complete, exit environment**:  
```bash
source deactivate ribosemap_env
```
