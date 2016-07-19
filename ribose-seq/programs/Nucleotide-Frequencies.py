#!/usr/bin/env python

#Author: Alli Gombolay
#This program calculates the ribonucleotide frequencies located at 5' position of input BED file

#Import Python modules
import sys
import os

#Module to create tables
from tabulate import tabulate

#Modules to create and read Excel files
import xlwt
import xlrd

#Module to parse command-line arguments
import argparse

#Use argparse function to create the "help" command-line option ([-h])
parser = argparse.ArgumentParser()
parser.add_argument('txt file')
parser.add_argument('subset')
parser.add_argument('location')
parser.add_argument("Location of user's local Ribose-seq directory")
args = parser.parse_args()

#Open input txt file and assign it to an object ("r": read file)
txt = open(sys.argv[1], "r")
subset = sys.argv[2]
location = sys.argv[3]
directory = sys.argv[4]

#Obtain name of input txt file excluding file extension
filename=os.path.splitext(os.path.basename(sys.argv[1]))[0]

#Obtain directory path of txt file
directory1=os.path.dirname(sys.argv[1])

#Specify directory path of output files
path="/".join(directory1.split('/'))
folder="/tables/"

#CALCULATE 5' NUCLEOTIDE FREQUENCIES

#Set the values of the base counts of nucleotide numbers to 0
A=0;
C=0;
G=0;
T=0;

for line in txt:
	if "A" in line:
		A+=1
	elif "C" in line:
		C+=1
	elif "G" in line:
		G+=1
	elif "T" in line:
		T+=1

#Calculate total number of nucleotides
total = (A+C+G+T)

#Calculate raw frequency of each nucleotide
A_frequency = float(A)/total
C_frequency = float(C)/total
G_frequency = float(G)/total
T_frequency = float(T)/total

#READ EXCEL FILE

background_frequencies = "%s/ribose-seq/results/Background-Nucleotide-Frequencies/%s.Nucleotide.Frequencies.xls" % (directory, subset)
workbook1 = xlrd.open_workbook(background_frequencies)
sheet1 = workbook1.sheet_by_index(0)

A_background = sheet1.cell_value(rowx=1, colx=2)
C_background = sheet1.cell_value(rowx=2, colx=2)
G_background = sheet1.cell_value(rowx=3, colx=2)
T_background = sheet1.cell_value(rowx=4, colx=2)

#Calculate normalized frequency of each nucleotide
A_normalized = A_frequency/A_background
C_normalized = C_frequency/C_background
G_normalized = G_frequency/G_background
T_normalized = T_frequency/T_background

#CREATE TABLE

#Create table of data with "tabulate" Python module
table = [["A", A, A_frequency, A_normalized, total], ["C", C, C_frequency, C_normalized, ""], ["G", G, G_frequency, G_normalized, ""], ["T", T, T_frequency, T_normalized, ""]]

#CREATE EXCEL FILE

#Create table of data with "xlwt" Python module
workbook2 = xlwt.Workbook()
sheet = workbook2.add_sheet("Sheet1")

decimal_style = xlwt.XFStyle()
decimal_style.num_format_str = '0.0000'

sheet.write(0, 0, "Ribonucleotide")
sheet.write(0, 1, "Number")
sheet.write(0, 2, "Raw Frequency")
sheet.write(0, 3, "Normalized Frequency")
sheet.write(0, 4, "Total")

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

sheet.write(1, 3, A_normalized, decimal_style)
sheet.write(2, 3, C_normalized, decimal_style)
sheet.write(3, 3, G_normalized, decimal_style)
sheet.write(4, 3, T_normalized, decimal_style)

sheet.write(1, 4, total)

#NAME OUTPUT FILES

#Specify name of output file based on input filename
output1=path+folder+filename+str('.Nucleotide.Frequencies.txt')
output2=path+folder+filename+str('.Nucleotide.Frequencies.xls')

#Redirect output to .txt file
sys.stdout=open(output1, "w")

#Specify header names and table style
print tabulate(table, headers=["Nucleotide", "Number", "Raw Frequency", "Normalized Frequency", "Total"], tablefmt="plain", floatfmt=".4f")

#Save table to .xls file
workbook2.save(output2)
