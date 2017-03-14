#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program determines outliers in the data using quartile analysis.

#Argument 1 = dataset file
dataset <- commandArgs(trailingOnly = TRUE)[1]

#Assign rNMP counts to variable
rNMPs <- read.csv(file=dataset, sep="\t")[ ,c('rNMPs')]

#First quartile
Q1 <- quantile(rNMPs, c(.25)) 

#Third quartile
Q3 <- quantile(rNMPs, c(.75)) 

#Interquartile range
IQR <- Q3-Q1

#Lower fence for outliers
lower.fence <- Q1-(3*IQR)

#Upper fence for outliers
upper.fence <- Q3+(3*IQR)

#Number of possible outliers in dataset
length <- length(rNMPs[which(rNMPs > upper.fence)])

#Print windows that are possible hotspots
read.csv(file=dataset, sep = "\t")[((length(rNMPs)-length+1):length(rNMPs)),]
