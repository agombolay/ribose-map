##Software Requirements:  
* [Python] (https://www.python.org/downloads/)  
 * [umitools] (https://github.com/brwnj/umitools): Trim UMIs and remove duplicate reads
 * [MACS2] (https://pypi.python.org/pypi/MACS2): Call peaks from alignment results

* [Bowtie] (https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/): Align sequencing reads to reference genome
 * [Manual page for information on commands] (http://bowtie-bio.sourceforge.net/manual.shtml)
 * [How to test if Bowtie index is properly installed] (http://bowtie-bio.sourceforge.net/tutorial.shtml)

* [SAMtools] (http://www.htslib.org/download/): Convert aligned reads files to BAM format
 * [Manual page for information on commands] (http://www.htslib.org/doc/samtools.html)
 * [Download SAMtools] (http://www.htslib.org/download/)

* [bedToBigBed] (http://hgdownload.cse.ucsc.edu/admin/exe/)

* [bedGraphToBigWig] (http://hgdownload.cse.ucsc.edu/admin/exe/): Convert files from bedGraph file format to bigWig format
  * bigWig file format is used to display dense, continuous data in the UCSC Genome Browser as a graph [(more)] (https://genome.ucsc.edu/goldenpath/help/bigWig.html)

* [bedtools]  (http://bedtools.readthedocs.org/en/latest/content/installation.html): Coverage tool calculates genome coverage based on alignment results
 * bedtools coverage tool calculates the depth and breadth of genome coverage from alignment results [(more)]
(http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html)

* [R]  (https://www.r-project.org/)
