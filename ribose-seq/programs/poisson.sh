#!/usr/bin/env bash

#Author: Alli Gombolay
#This program bins rNMPs into 2.5kb windows in reference genome

#Separate genome into 2.5 kb windows
bedtools makewindows -g /projects/home/agombolay3/data/repository/Ribose-seq-Project/ribose-seq/reference/sacCer2.bed -w 2500 > sacCer2.windows.bed

#Sort BED file of rNMP coordinates (sort first and second columns of file)
sort -k1,1 -k2,2n FS15.trimmed.v1.rNMP-coordinates.0-based.nuclear.txt > FS15.sorted.bed

#Determine intersect regions of the two BED files and then bin data
bedtools intersect -a sacCer2.windows.bed -b FS15.sorted.bed -c -sorted -nonamecheck > binned.data.txt
