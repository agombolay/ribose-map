#Set-up and use the Ribose-seq Analysis Pipeline
Author: Alli Gombolay  
Date: June 9, 2016

###1. [Install required software] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/Install-Software.md)

###2. [Install required modules] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/Install-Modules.md)

###3. Download raw sequencing files
* Sequencing files may be downloaded from [Illumina] (https://basespace.illumina.com/home/index) or [NCBI] (http://www.ncbi.nlm.nih.gov/gds/)

###4. [Convert sequencing files from .sra to .fq] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/Set-Up-Data-Files.sh)
* Sequencing files must be converted to FASTQ file format if not already

###5. [Download reference genome files from UCSC] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/Set-Up-Reference-Files.sh)

###6. [Analyze data with Ribose-seq Analysis programs] (https://github.com/agombolay/Ribose-seq-Project/tree/master/ribose-seq/scripts)
* The programs are numbered according to the order in which they should be run
* To view usage statements, please run the command below for any given program
```
<script name> -h
```

###Optional:
* If analyzing *S. cerevisiae*, [download expression data from Jay Hesselberth's site] (http://amc-sandbox.ucdenver.edu/User13/outbox/2016/)
