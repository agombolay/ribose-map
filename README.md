![Logo](https://github.com/agombolay/Images/blob/master/Logo.png)
# A bioinformatics toolkit for mapping rNMPs in genomic DNA
**© 2017 Alli Gombolay, Fredrik Vannberg, and Francesca Storici**  
**School of Biological Sciences, Georgia Institute of Technology**

**If you use Ribose-Map, please use the following citation**:  
Gombolay, AL, FO Vannberg, and F Storici. Ribose-Map: a bioinformatics toolkit to map ribonucleotides embedded in genomic DNA. *Nucleic Acids Research*, Volume 47, Issue 1, 10 Jan 2019, Page e5, https://doi.org/10.1093/nar/gky874. 

## Modules
1. **Alignment** aligns reads to reference genome using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and de-depulicates/de-multiplexes reads if needed

2. **Coordinate** calculates single-nucleotide genomic coordinates of rNMPs based on aligned reads for any currently available rNMP sequencing technique: ribose-seq, emRiboSeq, RHII-HydEn-seq, Alk-HydEn-seq, and Pu-seq.  

3. **Composition**: Calculates percentage of r{A, C, G, U} normalized to corresponding percentages of reference genome

4. **Sequence**: Calculates frequencies of r{A, C, G, U/T} at rNMP sites and up/downstream from those sites  

5. **Distribution**: Creates bedgraph files of per-nucleotide rNMP coverage to be visualized in a genome browser and plots the per-nucleotide rNMP coverage for each genomic unit of reference genome and normalizes for sequencing depth

6. **Hotspot**: Calculates top 1% most abundant sites of rNMP incorporation and creates sequence logos using [ggseqlogo](https://github.com/omarwagih/ggseqlogo)

## Required Data
1. Single or paired-end rNMP-seq data (FASTQ file)
2. Nucleotide sequence of reference genome (FASTA and FAI files)
3. Chromosome sizes of reference genome (extension .chrom.sizes)
4. Bowtie2 indexes of reference genome (extension = .bt2)
5. [Configuration file](https://github.com/agombolay/ribose-map/blob/master/lib/sample.config)

**Note**: If you have a BED file of single-nucleotide genomic coordinates and want to input that file directly into the Sequence, Distribution, and Hotspot Modules, create a folder with the filepath below, save the file in this folder, and create a config file.

"$repository/results/$sample/coordinate$quality" ($variables should be those provided in config)

## Install Software

1. **Create conda software environment**:  
* If you do not already have conda installed, please visit [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
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
ribose-map/modules/ribosemap {alignment, coordinate, composition, sequence, distribution, hotspot} config
```

3. **Once the analysis is complete, exit environment**:  
```bash
source deactivate ribosemap_env
```
