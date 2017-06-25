<h1><p align="center">Biological Information</p></h1>

### [*Escherichia coli K12*](https://www.ncbi.nlm.nih.gov/genome/167)
* GC Content: 50.60%
* Genome Size: 5.17 Mb
* Number of Chromosomes: 1
* [Download Genome from Illumina](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

### [*Homo sapiens*](https://www.ncbi.nlm.nih.gov/genome/51)
* GC Content: 40.90%
* Genome Size: 2996.43 Mb
* Number of Chromosomes: 23
* [Download Genome from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

### [*Mus musculus*](https://www.ncbi.nlm.nih.gov/genome/?term=Mus+musculus+%28house+mouse%29)
* GC Content: 42.58%
* Genome Size: 2,641.99 Mb
* Number of Chromosomes: 21
* [Download Genome from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/)

### [*Saccharomyces cerevisiae*](https://www.ncbi.nlm.nih.gov/genome/15)
* GC Content: 38.42%
* Genome Size: 12.30 Mb
* Number of Chromosomes: 16
* [Download Genome from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/)
* [Download Genome from SGD](http://www.yeastgenome.org/)

### [*Schizosaccharomyces pombe*](https://www.ncbi.nlm.nih.gov/genome/14)
* GC Content: 36.04%
* Genome Size: 12.59 Mb
* Number of Chromosomes: 3
* [Download Genome from PomBase](https://www.pombase.org/downloads/genome-datasets)

### [*Thermococcus barophilus*]
[Genomic sequence]:(https://www.ncbi.nlm.nih.gov/nuccore/NC_014804.1)
[Plasmid sequence]:(https://www.ncbi.nlm.nih.gov/nuccore/NC_015471.1)
* GC Content: 41.75%
* Genome Size: 2.23 Mb
    * One circular chromosome (2,020,078 bp)
    * One circular plasmid (54,159 bp)
* [Download Genome from Illumina](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

### [*Haloferax volcanii*](https://www.ncbi.nlm.nih.gov/genome/167)
* GC Content: 50.60%
* Genome Size: 5.17 Mb
* Number of Chromosomes: 1
* [Download Genome from Illumina](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

```
#Convert 2bit file to FASTA
twoBitToFa sacCer2.2bit sacCer2.fa

#Index FASTA file
samtools faidx sacCer2.fa
```
