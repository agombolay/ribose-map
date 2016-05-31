#!/usr/bin/env bash

#Author: Alli Gombolay
#This program automatically runs input FASTQ.GZ files through the Ribose-seq analysis pipeline
#Adapted from Jay Hesselberth's code located at https://github.com/hesselberthlab/modmap/tree/snake

if [ "$1" == "-h" ]; then
        exit;
fi
