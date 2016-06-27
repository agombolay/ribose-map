#!/usr/bin/env python

#Author: Alli Gombolay
#This program calculates the nucleotide frequencies of an input FASTA file

#Import Python modules
import sys
import os
import argparse

from tabulate import tabulate

#Use argparse function to create the "help" command-line option ([-h])
parser = argparse.ArgumentParser()
parser.add_argument('FASTA file')
args = parser.parse_args()

#Open input FASTA file and assign it to an object ("r": read file)
fasta = open(sys.argv[1], "r")

#CALCULATE NUCLEOTIDE FREQUENCIES

#Set the values of the base counts of nucleotide numbers to 0
A=0;
C=0;
G=0;
T=0;

for line in fasta:
    
	#Skip lines in file that start with ">" symbol (i.e., >2micron)
	if ">" in line:
		continue

	for character in line:
		if character == "A":
			A+=1
		if character == "C":
			C+=1
		if character == "G":
			G+=1
		if character == "T":
			T+=1

#Calculate total number of nucleotides
total = (A+C+G+T)

#Calculate frequency of each nucleotide
A_frequency = float(A)/total
C_frequency = float(C)/total 
G_frequency = float(G)/total
T_frequency = float(T)/total

#CREATE TABLE

#Create table of data with "tabulate" Python module
table = [["A", A, A_frequency, total], ["C", C, C_frequency, ""], ["G", G, G_frequency, ""], ["T", T, T_frequency, ""]]

#NAME OUTPUT FILE

#Obtain name of FASTA file without .fa file extension
filename=os.path.splitext(sys.argv[1])[0]

#Specify name of output file based on input filename
output=filename+str('.Nucleotide_Frequencies.txt')

#Redirect output to .txt file
sys.stdout=open(output, "w")

#Specify header names and table style
print tabulate(table, headers=["Nucleotide", "Number", "Frequency", "Total"], tablefmt="simple")
