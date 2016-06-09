##Description of Software Requirements:  

* [**Python**] (https://www.python.org/downloads/)  

   * Python Modules:

      * Modules contained within the Python Standard Library:
         * [itertools] (https://docs.python.org/2/library/index.html): Creates iterators for looping
         * [operator] (https://docs.python.org/2/library/index.html): Performs standard mathimatical operations
         * [collections] (https://docs.python.org/2/library/index.html): Employs specialized container data types

      * Additional Required Modules that need to be installed:
         * [pandas] (https://pypi.python.org/pypi/pandas/0.18.1/): Data analysis tools
         * [pyfaidx] (https://pypi.python.org/pypi/pyfaidx): Manipulation of FASTA files
         * [toolshed] (https://pypi.python.org/pypi/toolshed): Manipulation of tabular data files 
         * [pybedtools] (https://pypi.python.org/pypi/pybedtools): Extension of BEDtools program

------

* [**twoBitToFa**] (https://genome.ucsc.edu/goldenpath/help/twoBit.html)
  * twoBitToFa is used to convert .2bit files to .fa format for further processing

* [**umitools**] (https://github.com/brwnj/umitools):
  * umitools is used to trim unique molecular identifiers (UMIs) and to remove duplicate reads

* [**MACS2**] (https://pypi.python.org/pypi/MACS2):
  * MACS2 is used to call peaks based on alignment results to determine locations of hotspots

* [**Bowtie**] (https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/):
   * Bowtie is used to align the raw sequencing reads to the reference genome of interest
   * Resources:
      * [Manual page for information on commands] (http://bowtie-bio.sourceforge.net/manual.shtml)
      * [How to test if Bowtie index is properly installed] (http://bowtie-bio.sourceforge.net/tutorial.shtml)

* [**SAMtools**] (http://www.htslib.org/download/):
  * SAMtools is used to convert the output SAM file of aligned reads to BAM file format
  * Resources:
    * [Manual page for information on commands] (http://www.htslib.org/doc/samtools.html)
    * [Filter reads based on forward/reverse strands] (https://www.biostars.org/p/14378/)
    * [Explanation of SAMtools flags (i.e., "4")] (http://broadinstitute.github.io/picard/explain-flags.html)

* [**BEDtools**]  (http://bedtools.readthedocs.org/en/latest/content/installation.html):
  * BEDtools coverage tool is used to calculate genome coverage based on alignment results [(more)] (http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html)
  * BEDtools is also used in conjunction with pybedtools to calculate nucleotide frequencies [(more)] (https://pythonhosted.org/pybedtools/)

* [**bedToBigBed**] (http://hgdownload.cse.ucsc.edu/admin/exe/):
  * bedToBigBed is used to convert files from BED file format to bigBed format
  * bigBed files are used to display annotation items in the UCSC Genome Browser as a graph [(more)] (http://genome.ucsc.edu/goldenPath/help/bigBed.html)
  * Compared to BED files, bigBed files are a faster way to display data in the UCSC Genome Browser

* [**bedGraphToBigWig**] (http://hgdownload.cse.ucsc.edu/admin/exe/):
  * bedGraphToBigWig is used to convert files from bedGraph file format to bigWig format
  * bigWig files are used to display dense, continuous data in the UCSC Genome Browser as a graph [(more)] (https://genome.ucsc.edu/goldenpath/help/bigWig.html)
  * Compared to BedGraph files, bigWig files are a faster way to display data in the UCSC Genome Browser

* [**R**] (https://www.r-project.org/):
 * R software is used to create figures summarizing results of analysis

* [**Trimmomatic**] (http://www.usadellab.org/cms/?page=trimmomatic)
 * Trimmomatic program is used to trim reads based on quality (optional)

------

####To install the required software, please click on the links below:
* [Instructions on how to install Python and non-standard modules] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/Install-Python-Modules.md)  
* [Instructions on how to install required all other required software] (https://github.com/agombolay/Ribose-seq-Project/blob/master/Set-Up-Instructions/Install-Required-Software.md)
