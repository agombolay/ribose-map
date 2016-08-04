#! /usr/bin/env bash

#Author: Alli Gombolay
#Date: August 4, 2016
#This program allows the user to output only certain chromosomes from input FASTA file

samtools faidx $fasta $chromosomes > $output
