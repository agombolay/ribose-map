#! /usr/bin/env python

#Author: Alli Gombolay
#This program calculates the nucleotide frequencies of an input FASTA file

import sys
import os
import argparse

#Use argparse function to create the "help" command-line option ([-h])
parser = argparse.ArgumentParser()
parser.add_argument('FASTA file')
args = parser.parse_args()

#Open input FASTA file and assign it to an object ("r": read file)
fasta = open(sys.argv[1], "r")

#Set the values of the base counts of nucleotide frequencies to 0
A=0;
C=0;
G=0;
T=0;

#If necessary, skip the first line (description line) of FASTA file
#fasta.readline()

for line in fasta:
    
	#Skip lines in file that start with ">" symbol (i.e., >2micron)
	if ">" in line:
		continue

	#Optional: Convert letters in file lowercase for consistency
	#line = line.lower()

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

#Obtain name of FASTA file without .fa file extension
filename=os.path.splitext(sys.argv[1])[0]

#Specify name of output file based on input filename
output=filename+str('.txt')

#Redirect output to .txt file (i.e., 2micron.txt)
sys.stdout=open(output, "w")

#Convert freqeuncies into strings so they can be printed
print "Number of A's: " + str(A)
print "Number of C's: " + str(C)
print "Number of G's: " + str(G)
print "Number of T's: " + str(T)

print "Total number of nucleotides: " + str(total)

#Calculate ratio of G's and C's to total base content
#0. = Convert number to a floating point decimal number
#gc = (G+C+0.) / (A+C+G+T+0.)

#Print GC content value
#print "GC content: " + str(gc)

print "Frequency of A's: " + str(float(A)/total)
print "Frequency of C's: " + str(float(C)/total)
print "Frequency of G's: " + str(float(G)/total)
print "Frequency of T's: " + str(float(T)/total)

