#Set-up and use the Ribose-seq Analysis Pipeline
Author: Alli Gombolay  
Date: June 9, 2016

###1. [Install required software] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/3_Software-Installation.md)

###2. [Install required modules] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/4_Modules-Installation.md)

###3.
a) Download raw sequencing files

b) [Convert sequencing files from .sra to .fq] (https://github.com/agombolay/Ribose-seq-Project/blob/master/ribose-seq/scripts/1_setUpRawData.sh)
* Sequencing files must be converted to FASTQ file format if not already

###4. [Download reference genome files] (https://github.com/agombolay/Ribose-seq-Project/blob/master/ribose-seq/scripts/1_setUp.sh)

###5. [Analyze data with Ribose-seq Analysis scripts] (https://github.com/agombolay/Ribose-seq-Project/tree/master/ribose-seq/scripts)
* The scripts are numbered according to the order in which they should be run

###Optional:
* If analyzing *S. cerevisiae*, [download expression data from Jay Hesselberth's site] (http://amc-sandbox.ucdenver.edu/User13/outbox/2016/)
