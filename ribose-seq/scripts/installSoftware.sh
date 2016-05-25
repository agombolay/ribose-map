#INSTALL MACS2 PROGRAM TO CALL ALIGNMENT PEAKS

#Download the compressed file from the Python MACS2 website
wget https://pypi.python.org/packages/9f/99/a8ac96b357f6b0a6f559fe0f5a81bcae12b98579551620ce07c5183aee2c/MACS2-2.1.1.20160309.tar.gz

#Uncompress the downloaded .tar.gz file
tar -zxvf MACS2-2.1.1.20160309.tar.gz

#Install the MACS2 program with Python version 2.7
#"--user": Allows user to avoid permission issues
python2.7 setup.py install --user

#Add the path to the MACS2 script to your $PATH via your .bashrc file
echo "export PATH="/projects/home/agombolay3/.local/bin:$PATH"" >> ~/.bashrc

#INSTALL UMITOOLS PROGRAM TO REMOVE PCR DUPLICATES

#Install the pysam program
#Pip must already be installed
pip install pysam

#Install the editdist program
#Download the compressed file from the editdist website
wget https://py-editdist.googlecode.com/files/py-editdist-0.3.tar.gz

#Uncompress the downloaded .tar.gz file
tar xzf py-editdist-0.3.tar.gz

#Change working directory to editdist folder
cd py-editdist-0.3/

#Install the editdist program with Python version 2.7
#"--user": Allows user to avoid permission issues
python2.7 setup.py install --user

#Download the compressed file from the umitools website
wget -O umitools-master.zip https://github.com/brwnj/umitools/archive/master.zip

#Uncompress the downloaded .tar.gz file
unzip umitools-master.zip

#Change working directory to umitools folder
cd umitools-master

#Install the umitools program with Python version 2.7
#"--user": Allows user to avoid permission issues
python2.7 setup.py install --user

#Add the path to the umitools script to your $PATH via your .bashrc file
echo "export PATH="/projects/home/agombolay3/.local/bin:$PATH"" >> ~/.bashrc
