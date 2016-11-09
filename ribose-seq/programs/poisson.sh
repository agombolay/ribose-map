#!/usr/bin/env bash

#Author: Alli Gombolay
#This program bins rNMPs into 2.5kb windows in reference genome

#Separate genome into 2.5 kb windows
bedtools makewindows -g chrom.sizes -w 2500 > genome.windows.bed

#Sort BED file of rNMP coordinates (sort first and second columns of file)
sort -k1,1 -k2,2n FS15.sorted.bed > temporary && mv temporary FS15.sorted.bed

#Determine intersect regions of the two BED files and then bin data
bedtools intersect -a genome.windows.bed -b FS15.sorted.bed -c -sorted
