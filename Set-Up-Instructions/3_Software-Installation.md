##Instructions to install required software

####SRA Toolkit
1. Download the SRA Toolkit program from NCBI's website
```
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.6.3/sratoolkit.2.6.3-centos_linux64.tar.gz
```

####twoBitToFa
1. Download the twoBitToFa program from UCSC's application site  
(Reference: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
```
&nbsp; 
2. Change file permissions to make the program executable
```
chmod +x twoBitToFa
```

####Bowtie
1. Download the .zip file from the bowtie website
```
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip
```
&nbsp; 
2. Uncompress the downloaded .zip file
```
unzip bowtie-1.1.2-linux-x86_64.zip
```
&nbsp;
3. Add the path to the umitools script to your $PATH via the .bashrc configuration file
```
echo "export PATH="/projects/home/agombolay3/data/bin/bowtie-1.1.2:$PATH"" >> ~/.bashrc
```

####SAMtools
1. Download the .zip file from the SAMtools website  
(Reference: http://www.htslib.org/download/)
```
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
```
&nbsp;
2. Uncompress the downloaded .tar.bz2 file
```
tar -jxvf samtools-1.3.1.tar.bz2
```
&nbsp;
3. Change current directory to samtools folder
```
cd samtools-1.3.1/
```
&nbsp;
4. Execute commands outlined in the Makefile
```
make
```
&nbsp;
5. Install the SAMtools program to specified location
```
make prefix=/projects/home/agombolay3/data/bin install
```
&nbsp;
6. Add the path to the umitools script to your $PATH via the .bashrc configuration file
```
echo "export PATH="/projects/home/agombolay3/data/bin/bin:$PATH"" >> ~/.bashrc
```

####BEDtools
1. Download the .tar.gz file from the BEDtools website 
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
```
&nbsp;
2. Uncompress the downloaded .tar.gz file
```
tar -zxvf bedtools-2.25.0.tar.gz
```
&nbsp;
3. Change current directory to bedtools folder
```
cd bedtools2
```
&nbsp;
4. Execute commands outlined in the Makefile
```
make
```
&nbsp;
5. Add the path to the umitools script to your $PATH via the .bashrc configuration file
```
echo "export PATH="/projects/home/agombolay3/data/bin/bedtools2/bin:$PATH"" >> ~/.bashrc
```

####umitools
(Note: Instructions adapted from https://github.com/brwnj/umitools)

####Part a.) Pysam Program

1. Install the pysam program  
* Pip must already be installed
```
pip install pysam
```

####Part b.) Editdist Program

1. Download the .tar.gz file from the editdist website  
(Reference: http://www.mindrot.org/projects/py-editdist/)
```
wget https://py-editdist.googlecode.com/files/py-editdist-0.3.tar.gz
```
&nbsp; 
2. Uncompress the downloaded .tar.gz file
```
tar xzf py-editdist-0.3.tar.gz
```
&nbsp; 
3. Change current directory to editdist folder
```
cd py-editdist-0.3/
```
&nbsp; 
4. Install the editdist program with Python version 2.7
* "--user" allows bypasses permission issues on server
```
python2.7 setup.py install --user
```

####Part c.) Umitools Program

1. Download the .zip file from the umitools website  
(Reference: https://github.com/brwnj/umitools)
```
wget -O umitools-master.zip https://github.com/brwnj/umitools/archive/master.zip
```
&nbsp; 
2. Uncompress the downloaded .zip file
```
unzip umitools-master.zip
```
&nbsp; 
3. Change current directory to umitools folder
```
cd umitools-master
```
&nbsp; 
4. Install the umitools program with Python version 2.7 
* "--user" allows bypasses permission issues on server
```
python2.7 setup.py install --user
```
&nbsp; 
5. Add the path to the umitools script to your $PATH via the .bashrc configuration file
```
echo "export PATH="/projects/home/agombolay3/.local/bin:$PATH"" >> ~/.bashrc
```

####MACS2

1. Download the .tar.gz file from the MACS2 website  
(Reference: https://pypi.python.org/pypi/MACS2)
```
wget https://pypi.python.org/packages/9f/99/a8ac96b357f6b0a6f559fe0f5a81bcae12b98579551620ce07c5183aee2c/MACS2-2.1.1.20160309.tar.gz
```
&nbsp; 
2. Uncompress the downloaded .tar.gz file
```
tar -zxvf MACS2-2.1.1.20160309.tar.gz
```
&nbsp; 
3. Install the MACS2 program with Python version 2.7  
* "--user" allows bypasses permission issues on server
```
python2.7 setup.py install --user
```
&nbsp; 
4. Add the path to the MACS2 script to your $PATH via the .bashrc configuration file
```
echo "export PATH="/projects/home/agombolay3/.local/bin:$PATH"" >> ~/.bashrc
```

###bedToBigBed
1. Download the program from the UCSC Genome Browswer website  
(Reference: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
```
&nbsp;
2. Change the file permissions to make the file executable
```
chmod +x bedToBigBed
```
&nbsp;
3. Add the path to the bedToBigBed script to your $PATH via the .bashrc configuration file
```
echo "export PATH="/projects/home/agombolay3/data/bin:$PATH"" >> ~/.bashrc
```

###bedgraphToBigWig
1. Download the program from the UCSC Genome Browswer website  
(Reference: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedgraphToBigWig
```
&nbsp;
2. Change the file permissions to make the file executable
```
chmod +x bedgraphToBigWig
```
&nbsp;
3. Add the path to the bedToBigBed script to your $PATH via the .bashrc configuration file
```
echo "export PATH="/projects/home/agombolay3/data/bin:$PATH"" >> ~/.bashrc
```
