#!/usr/bin/env python

#Author: Alli Gombolay
#This program calculates the nucleotide frequencies of an input FASTA file

#Import Python modules
import sys
import os

#Module to create tables
from tabulate import tabulate

#Module	to create Excel	files
import xlwt

#Module	to parse command-line arguments
import argparse

#Use argparse function to create command-line options
class SmartFormatter(argparse.HelpFormatter):

	def _split_lines(self, text, width):
		if text.startswith('R|'):
			return text[2:].splitlines()
		return argparse.HelpFormatter._split_lines(self, text, width)

from argparse import ArgumentParser

parser = ArgumentParser(description='This program calculates background nucleotide Frequencies; Requires Python2.7+', formatter_class=SmartFormatter)

parser.add_argument('-i', choices=['FASTA'], help="R|Filepath of FASTA file (/path/to/reference.fasta)")

parser.parse_args()

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

#CREATE EXCEL FILE

#Create table of data with "xlwt" Python module
workbook = xlwt.Workbook()
sheet = workbook.add_sheet("Sheet1")

decimal_style = xlwt.XFStyle()
decimal_style.num_format_str = '0.0000'

sheet.write(0, 0, "Nucleotide")
sheet.write(0, 1, "Number")
sheet.write(0, 2, "Frequency")
sheet.write(0, 3, "Total")

sheet.write(1, 0, "A")
sheet.write(2, 0, "C")
sheet.write(3, 0, "G")
sheet.write(4, 0, "T")

sheet.write(1, 1, A)
sheet.write(2, 1, C)
sheet.write(3, 1, G)
sheet.write(4, 1, T)

sheet.write(1, 2, A_frequency, decimal_style)
sheet.write(2, 2, C_frequency, decimal_style)
sheet.write(3, 2, G_frequency, decimal_style)
sheet.write(4, 2, T_frequency, decimal_style)

sheet.write(1, 3, total)

#NAME OUTPUT FILE

#Obtain name of input FASTA file excluding file extension
filename=os.path.splitext(os.path.basename(sys.argv[1]))[0]

#Obtain directory path of FASTA file
directory=os.path.dirname(sys.argv[1])

#Extract only a part of directory path
path="/".join(directory.split('/')[:-1])

#Specify directory path of output files
folder="/results/Background-Nucleotide-Frequencies/"

#Specify name of output file based on input filename
output1=path+folder+filename+str('.Nucleotide.Frequencies.txt')
output2=path+folder+filename+str('.Nucleotide.Frequencies.xls')

#Redirect output to .txt file
sys.stdout=open(output1, "w")

#Specify header names and table style
print tabulate(table, headers=["Nucleotide", "Number", "Frequency", "Total"], tablefmt="simple", floatfmt=".4f")

#Save table to .xls file
workbook.save(output2)
