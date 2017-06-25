<h1><p align="center">Biological Information</p></h1>

### [*Escherichia coli K12*](https://www.ncbi.nlm.nih.gov/genome/167)
* GC content: 50.60%
* Genome size: 5.17 Mb
* Number of chromosomes: 1
* [Download genome from Illumina](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

### [*Homo sapiens*](https://www.ncbi.nlm.nih.gov/genome/51)
* GC content: 40.90%
* Genome size: 2996.43 Mb
* Number of chromosomes: 23
* [Download genome from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

### [*Mus musculus*](https://www.ncbi.nlm.nih.gov/genome/?term=Mus+musculus+%28house+mouse%29)
* GC content: 42.58%
* Genome size: 2,641.99 Mb
* Number of chromosomes: 21
* [Download genome from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/)

### [*Saccharomyces cerevisiae*](https://www.ncbi.nlm.nih.gov/genome/15)
* GC content: 38.42%
* Genome size: 12.30 Mb
* Number of chromosomes: 16
* [Download genome from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/)
* [Download genome from SGD](http://www.yeastgenome.org/)

### [*Schizosaccharomyces pombe*](https://www.ncbi.nlm.nih.gov/genome/14)
* GC content: 36.04%
* Genome size: 12.59 Mb
* Number of chromosomes: 3
* [Download genome from PomBase](https://www.pombase.org/downloads/genome-datasets)

### [*Thermococcus barophilus*](https://www.ncbi.nlm.nih.gov/genome/1564)
* GC content: 41.75%
* Genome size: 2.23 Mb
     * One circular plasmid (54,159 bp)
     * One circular chromosome (2,020,078 bp)
* [Download genomic sequence from NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_014804.1)
* [Download plasmid sequence from NCBI:](https://www.ncbi.nlm.nih.gov/nuccore/NC_015471.1)

### [*Haloferax volcanii*](https://www.ncbi.nlm.nih.gov/genome/167)
* GC content: 50.60%
* Genome size: 5.17 Mb
* Number of chromosomes: 1
* [Download genome from Illumina](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

```
#Convert 2bit file to FASTA
twoBitToFa sacCer2.2bit sacCer2.fa

#Index FASTA file
samtools faidx sacCer2.fa
```
