#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/extract.fq
