#!/usr/bin/env bash

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

umi_tools extract -v 0 -I $fastq1 -p $UMI -S $output/extract.fq

grep -B 1 -A 2 ^$barcode $output/extract.fq | sed '/^--$/d' \
| awk 'NR%2 == 0 {sub(/^.{'${#barcode}'}/,"")} {print}' > $output/filter.fq
