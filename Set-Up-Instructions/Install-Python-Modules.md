#Instructions to install Python

(Note: Instructions adapted from [here]) (https://mail.python.org/pipermail/tutor/2002-March/012903.html)

```
wget https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz
```

```
tar zxf Python-2.7.11.tgz
```

```
cd Python-2.2
```

```
./configure
```

```
make
```

------

#Instructions to install modules

###References
* [Installing Python modules] (https://docs.python.org/2/install/)
* [Using "--user" flag to install] (http://stackoverflow.com/questions/18199853/error-could-not-create-library-python-2-7-site-packages-xlrd-permission-den)

###[Toolshed] (https://pypi.python.org/pypi/toolshed)

```
wget https://pypi.python.org/packages/2e/97/f9a1a64d70e47d1c0d84a1a6671ba7828b7f99df09ab83103b8cce838406/toolshed-0.4.5.tar.gz
```
```
tar -zxvf toolshed-0.4.5.tar.gz 
```
```
cd toolshed-0.4.5
```
```
python setup.py install --user
```

###[Pyfaidx] (https://pypi.python.org/pypi/pyfaidx)

```
wget https://pypi.python.org/packages/9f/e6/e745e065c291b25e3fd025360edfe6289c2ee35ffa9375e85708612ca820/pyfaidx-0.4.7.1.tar.gz
```
```
tar -zxvf pyfaidx-0.4.7.1.tar.gz 
```
```
cd pyfaidx-0.4.7.1
```
```
python setup.py install --user
```

###[Pybedtools] (https://pypi.python.org/pypi/pybedtools)

```
wget https://pypi.python.org/packages/c2/a6/29fd41e86b38e7ca39eb7818c422e927534f2bd324c34cfd411cfca203b0/pybedtools-0.7.7.tar.gz
```
```
tar -zxvf pybedtools-0.7.7.tar.gz 
```
```
cd pybedtools-0.7.7
```
```
python setup.py install --user
```

###[Pandas] (https://pypi.python.org/pypi/pandas)

```
wget https://pypi.python.org/packages/11/09/e66eb844daba8680ddff26335d5b4fead77f60f957678243549a8dd4830d/pandas-0.18.1.tar.gz
```
```
tar -zxvf pandas-0.18.1.tar.gz
```
```
cd pandas-0.18.1
```
```
python setup.py install --user
```
