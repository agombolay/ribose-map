#!/usr/bin/env python

#Author: Alli Gombolay
#This program calculates the nucleotide frequencies of the 5' position of an input FASTA file

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

#Obtain name of input FASTA file excluding file extension
filename=os.path.splitext(os.path.basename(sys.argv[1]))[0]

#Obtain directory path of FASTA file
directory=os.path.dirname(sys.argv[1])

#Specify directory path of output files
path="/".join(directory.split('/')[:-1])
folder="/nucleotideFrequencies/"

#Specify name of output file based on input filename
list=path+folder+filename+str('.5-Prime-Nucleotides-List.txt')

#Create file to where list of 5' nucleotides will saved
standard_output = sys.stdout
list1 = open(list, 'w')
sys.stdout = list1

#CALCULATE 5' NUCLEOTIDE FREQUENCIES

#Skip header lines and print only 5' nucleotides (line[:1])
for line in fasta:

        if ">" in line:
                continue
	
	print line[:1]

#Redirect standard output to file of list of 5' nucleotides
sys.stdout = standard_output

#Close the file
list1.close()

#Open list and assign it to a different object ("r": read file)
list2 = open(list, 'r')

#Set the values of the base counts of 5' nucleotide numbers to 0
A=0;
C=0;
G=0;
T=0;

for line in list2:

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
table = [["A", T, T_frequency, total], ["C", G, G_frequency, ""], ["G", C, C_frequency, ""], ["U", A, A_frequency, ""]]

#NAME OUTPUT FILE

#Specify name of output file based on input filename
output=path+folder+filename+str('.Ribonucleotide.Frequencies.txt')

#Redirect output to .txt file
sys.stdout=open(output, "w")

#Specify header names and table style
print tabulate(table, headers=["Nucleotide", "Number", "Frequency", "Total"], tablefmt="simple")
