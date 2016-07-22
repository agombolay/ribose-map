#!/usr/bin/env python

#Author: Alli Gombolay
#This program calculates the ribonucleotide frequencies located at 3' position of input BED file

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
parser.add_argument('BED file')
parser.add_argument('Reference')
parser.add_argument("Location of user's local Ribose-seq directory")
args = parser.parse_args()

#Open input BED file and assign it to an object ("r": read file)
bed = open(sys.argv[1], "r")
reference = sys.argv[2]
directory1 = sys.argv[3]

#Obtain name of input BED file excluding file extension
filename=os.path.splitext(os.path.basename(sys.argv[1]))[0][:-12]

#Obtain directory path of BED file
directory2=os.path.dirname(sys.argv[1])

#Specify directory path of output files
path="/".join(directory2.split('/')[:-1])
folder="/Ribonucleotides/tables/"

#Specify name of output file based on input filename
list=path+folder+filename+str('.')+reference+str('.Ribonucleotides-List.txt')

#Create file to where list of ribonucleotides will saved
standard_output = sys.stdout
list1 = open(list, 'w')
sys.stdout = list1

#CALCULATE 5' NUCLEOTIDE FREQUENCIES

#Print only ribonucleotides (3' end of read (end for + strand and start for - strand))
for line in bed:
	
	if (reference == "sacCer2" and "+" in line):
		print(line.split()[4])[-1],
		print(line.split()[5])

	elif (reference == "sacCer2" and "-" in line):
		print(line.split()[4])[:1],
		print(line.split()[5])

	elif (reference == "nuclear" and "chrM" not in line and "+" in line):
		print(line.split()[4])[-1],
		print(line.split()[5])

	elif (reference == "nuclear" and "chrM" not in line and "-" in line):
                print(line.split()[4])[:1],
		print(line.split()[5])

	elif (reference == "chrM" and "chrM" in line and "+" in line):
		print(line.split()[4])[-1],
		print(line.split()[5])

	elif (reference == "chrM" and "chrM" in line and "-" in line):
		print(line.split()[4])[:1],
		print(line.split()[5])

	elif (reference == "2micron" and "2micron" in line  and "+" in line):
		print(line.split()[4])[-1],
		print(line.split()[5])

	elif (reference == "2micron" and "2micron" in line and "-" in line):
                print(line.split()[4])[:1],
		print(line.split()[5])

#Redirect standard output to file of list of 3' nucleotides
sys.stdout = standard_output

#Close the file
list1.close()

#Open list and assign it to a different object ("r": read file)
list2 = open(list, 'r')

#Set the values of the base counts of 3' nucleotide numbers to 0
A=0;
C=0;
G=0;
U=0;

for line in list2:

	#for character in line:

	if ("A" in line and "+" in line):
		A+=1
	elif ("A" in line and "-" in line):
		U+=1
	elif ("C" in line and "+" in line):
		C+=1
	elif ("C" in line and "-" in line):
		G+=1
	elif ("G" in line and "+" in line):
		G+=1
	elif ("G" in line and "-" in line):
		C+=1
	elif ("T"  in line and "+" in line):
		U+=1
	elif ("T" in line and "-" in line):
		A+=1

#Calculate total number of nucleotides
total = (A+C+G+U)

#Calculate raw frequency of each nucleotide
A_frequency = float(A)/total
C_frequency = float(C)/total
G_frequency = float(G)/total
U_frequency = float(U)/total

#READ EXCEL FILE

background_frequencies = "%s/ribose-seq/results/Background-Nucleotide-Frequencies/%s.Nucleotide.Frequencies.xls" % (directory1, reference)
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
U_normalized = U_frequency/T_background

#CREATE TABLE

#Create table of data with "tabulate" Python module
table = [["A", A, A_frequency, A_normalized, total], ["C", C, C_frequency, C_normalized, ""], ["G", G, G_frequency, G_normalized, ""], ["U", U, U_frequency, U_normalized, ""]]

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
sheet.write(4, 0, "U")

sheet.write(1, 1, A)
sheet.write(2, 1, C)
sheet.write(3, 1, G)
sheet.write(4, 1, U)

sheet.write(1, 2, A_frequency, decimal_style)
sheet.write(2, 2, C_frequency, decimal_style)
sheet.write(3, 2, G_frequency, decimal_style)
sheet.write(4, 2, U_frequency, decimal_style)

sheet.write(1, 3, A_normalized, decimal_style)
sheet.write(2, 3, C_normalized, decimal_style)
sheet.write(3, 3, G_normalized, decimal_style)
sheet.write(4, 3, U_normalized, decimal_style)

sheet.write(1, 4, total)

#NAME OUTPUT FILES

#Specify name of output file based on input filename
output1=path+folder+filename+str('.')+reference+str('.Ribonucleotide.Frequencies.txt')
output2=path+folder+filename+str('.')+reference+str('.Ribonucleotide.Frequencies.xls')

#Redirect output to .txt file
sys.stdout=open(output1, "w")

#Specify header names and table style
print tabulate(table, headers=["Ribonucleotide", "Number", "Raw Frequency", "Normalized Frequency", "Total"], tablefmt="simple", floatfmt=".4f")

#Save table to .xls file
workbook2.save(output2)
