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
make_option(c("-c", "--forward"), help="Number of observed ribonucleotides on forward strand"),
make_option(c("-c", "--reverse"), help="Number of observed ribonucleotides on reverse strand"),
make_option(c("-r", "--reference"), help="Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)"),
make_option(c("-d", "--directory"), help="Local user directory (e.g., /projects/home/agombolay3/data/repository)")
)

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

for(strand in c("Forward", "Reverse")) {
	
	#Specify output directory and file
	output <- file.path(opt$directory, "Results", opt$reference, opt$sample, "Random")
	file <- file.path(output, paste(opt$sample, "-", "Random", ".", strand, ".txt", sep=""))
	
	data = c()
	for (i in 1:opt$cells) {
		
		if (strand == "Forward") {

		data = sample(1:opt$positions, opt$forward, replace=TRUE)
		vector <- append(data,sample(1:opt$positions, opt$forward, replace=FALSE))
		
		}
		
		else if (strand == "Reverse") {

		data = sample(1:opt$positions, opt$reverse, replace=TRUE)
		vector <- append(data,sample(1:opt$positions, opt$reverse, replace=FALSE))
		
		}
	}
}

hx=hist(vector,breaks=seq(1,opt$positions,l=opt$positions+1),plot=FALSE)
plot(hx$counts,xlab='Position on chromosome',ylab='Frequency of rNMPs')
