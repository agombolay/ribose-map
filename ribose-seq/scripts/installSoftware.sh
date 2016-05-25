#INSTALL MACS2 PROGRAM TO CALL ALIGNMENT PEAKS

#Download file from MACS2 website
wget https://pypi.python.org/packages/9f/99/a8ac96b357f6b0a6f559fe0f5a81bcae12b98579551620ce07c5183aee2c/MACS2-2.1.1.20160309.tar.gz

#Uncompress downloaded .tar.gz file
tar -zxvf MACS2-2.1.1.20160309.tar.gz

#Install MACS2 program with Python version 2.7
#"--user": Allows user to avoid permission issues
python2.7 setup.py install --user

#Add the path to the MACS2 script to your $PATH via your .bashrc file
echo "export PATH="/projects/home/agombolay3/.local/bin:$PATH"" >> ~/.bashrc
