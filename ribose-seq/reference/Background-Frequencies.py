#!/usr/bin/env python

#Author: Alli Gombolay
#This program calculates the nucleotide frequencies of an input FASTA file

#Import Python modules
import sys
import os
import argparse
#import tablib
from tabulate import tabulate

#Use argparse function to create the "help" command-line option ([-h])
parser = argparse.ArgumentParser()
parser.add_argument('FASTA file')
args = parser.parse_args()

#Open input FASTA file and assign it to an object ("r": read file)
fasta = open(sys.argv[1], "r")

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

#CREATE EXCEL FILE

#data = tablib.Dataset()

#Add data values to Excel file
#data.append_col(['A', 'C', 'G', 'T'], header='Nucleotide')
#data.append_col([A, C, G, T], header='Number of Occurrences')
#data.append_col([A_frequency, C_frequency, G_frequency, T_frequency], header='Frequency')
#data.append_col([total, '     ', '     ', '     '], header='Total Number of Nucleotides')

#Add column headers to Excel file
#data.headers = ['Nucleotide', 'Number', 'Frequency', 'Total Number of Nucleotides']

#NAME EXCEL FILE

#Obtain name of FASTA file without .fa file extension
#filename=os.path.splitext(sys.argv[1])[0]

#Specify name of output file based on input filename
#output=filename+str('.Nucleotide_Frequencies.xls')

#Redirect output to .xls file
#sys.stdout=open(output, "w")

#print data.xls

table = [["A", A, A_frequency, total], ["C", C, C_frequency, ""], ["G", G, G_frequency, ""], ["T", T, T_frequency, ""]]

#Obtain name of FASTA file without .fa file extension
filename=os.path.splitext(sys.argv[1])[0]

#Specify name of output file based on input filename
output=filename+str('.Nucleotide_Frequencies.txt')

sys.stdout=open(output, "w")

print tabulate(table, headers=["Nucleotide", "Number", "Frequency", "Total"], tablefmt="simple")
