#!/usr/bin/env Rscript

#Copyright 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots rNMP and flanking nucleotide frequencies (100 bp upstream and downstream)
#rNMP frequencies = position 0, upstream frequencies = -100 -> -1, and downstream = +1 -> +100

#############################################################################################################################
#Load libraries
library(optparse)
library(ggplot2)

#Command line options
option_list <- list(make_option(c("-i", "--input"), help="Input sample name"),
make_option(c("-t", "--title"), help="Title will be same for all plots; blank = no title"),
make_option(c("-r", "--reference"), help="Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)"),
make_option(c("-d", "--directory"), help="Local user directory (e.g., /projects/home/agombolay3/data/repository)"))

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

#############################################################################################################################
for(i in c("all", "mito", "nucleus")) {

    #Specify output directory and file
    output <- file.path(opt$directory, "Ribose-Map", "Results", opt$reference, opt$input, "Frequencies")
    file <- file.path(output, paste(opt$input, "-", "Frequencies", ".", opt$reference, ".", i, ".txt", sep=""))

        #Plot only if files exist
            if (file.exists(file)) {
		#Plot regular and zoomed datasets
                for(j in c("regular", "zoomed")) {

#############################################################################################################################
                    #Specify datasets to be used for each round of loop
                    if (j == "regular") {data = read.table(file, sep="\t", header=TRUE)}
                    if (j == "zoomed") {data = read.table(file, sep="\t", header=TRUE)[86:116,]}
    
                    #Define variables to store nucleotide positions and frequency values
                    position <- data$X; A1 <- data$A; C1 <- data$C; G1 <- data$G; T1 <- data$U.T

#############################################################################################################################
                    #Plot frequencies
                    myplot <- ggplot(data=data, aes(x=position)) +
    
                    #Plot data as scatterplot
                    geom_line(aes(y = A1, colour = "A")) + geom_point(aes(y = A1, colour = "A")) +
                    geom_line(aes(y = C1, colour = "C")) + geom_point(aes(y = C1, colour = "C")) +
                    geom_line(aes(y = G1, colour = "G")) + geom_point(aes(y = G1, colour = "G")) +
                    geom_line(aes(y = T1, colour = "U/T")) + geom_point(aes(y = T1, colour = "U/T")) +
    
                    #Add axes titles and plot title
                    xlab("Position") + ylab("Frequency") + ggtitle(opt$title) +
    
                    #Increase distance between axes and lables
                    theme(axis.title.y=element_text(margin=margin(0,20,0,0))) +
                    theme(axis.title.x=element_text(margin=margin(20,0,0,0))) +
    
                    #Remove and replace default background plot theme
                    theme_bw() + theme(panel.border=element_blank(), panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) +
    
                    #Specify font size for plot text and center title of plot
                    theme(text = element_text(size=14)) + theme(plot.title = element_text(hjust = 0.5)) +

                    #Specify color values for each of the four different nucleotides
                    scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="Nucleotide")

#############################################################################################################################
                    #Specify output path and save plot as PNG file
                    ggsave(filename=file.path(output, paste(opt$input, "-", j, ".", i, ".png", sep="")), plot=myplot)
                }
            }
}
