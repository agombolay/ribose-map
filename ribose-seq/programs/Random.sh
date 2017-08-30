#!/usr/bin/env Rscript

#© 2017 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots observed distribution of ribonucleotides vs random distribution

#Load library
library(optparse)

#Command line options
option_list <- list(
make_option(c("-c", "--cells"), help="Estimated number of cells"),
make_option(c("-s", "--sample"), help="Sample name(s) (e.g., FS1, FS2, FS3)"),
make_option(c("-r", "--reference"), help="Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)"),
make_option(c("-d", "--directory"), help="Local user directory (e.g., /projects/home/agombolay3/data/repository)")
)

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

data=c()
for (i in 1:cells) {
	data=sample(1:positions, ribos, replace=TRUE)
	vector <- append(data,sample(1:positions, ribos, replace=FALSE))
}

hx=hist(vector,breaks=seq(1,positions,l=positions+1),plot=FALSE)
plot(hx$counts)
