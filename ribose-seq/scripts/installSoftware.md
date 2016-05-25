###INSTALL MACS2 PROGRAM TO CALL ALIGNMENT PEAKS

1. Download the compressed file from the Python MACS2 website (https://pypi.python.org/pypi/MACS2)
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

###INSTALL UMITOOLS PROGRAM TO REMOVE PCR DUPLICATES
(Note: Instruction adapted from https://github.com/brwnj/umitools)

####Part a.) Pysam Program

Install the pysam program  
* Pip must already be installed
```
pip install pysam
```

####Part b. Editdist Program

**Install the editdist program**

1. Download the compressed file from the editdist website
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

**Install the umitools program**

1. Download the compressed file from the umitools website
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
5. Add the path to the umitools script to your $PATH via your .bashrc file
```
echo "export PATH="/projects/home/agombolay3/.local/bin:$PATH"" >> ~/.bashrc
```
