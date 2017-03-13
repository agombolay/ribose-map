#!/usr/bin/env Rscript

#Copyright 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots rNMP and flanking nucleotide frequencies (100 bp upstream and downstream)
#rNMP frequencies = position 0, upstream frequencies = -100 -> -1, and downstream = +1 -> +100

#Define three command line arguments
file <- commandArgs(trailingOnly = TRUE)[1]
title <- commandArgs(trailingOnly = TRUE)[2]
filename <- commandArgs(trailingOnly = TRUE)[3]

if (file.exists(file)) {
    #Define data file (tab delimited) to read
    data <- read.table(file, sep="\t", header=TRUE)

    #Define position values
    position <- data$X
    
    #Define frequency values
    frequencyA <- data$A; frequencyC <- data$C
    frequencyG <- data$G; frequencyT <- data$U.T

    #Load ggplot2
    library(ggplot2)

    #Plot frequencies
    myplot <- ggplot(data=data, aes(x=position)) +
        
        #Plot data as scatterplot   
        geom_line(aes(y = frequencyA, colour = "A")) + geom_point(aes(y = frequencyA, colour = "A")) +
        geom_line(aes(y = frequencyC, colour = "C")) + geom_point(aes(y = frequencyC, colour = "C")) +
        geom_line(aes(y = frequencyG, colour = "G")) + geom_point(aes(y = frequencyG, colour = "G")) +
        geom_line(aes(y = frequencyT, colour = "U/T")) + geom_point(aes(y = frequencyT, colour = "U/T")) +
    
        #Remove and replace default background plot theme
        theme_bw() + theme(panel.border=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) +
    
        #Specify color values for each of the four nucleotides
        scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="Nucleotide") +
    
        #Label axes, add title, define y-axis limits, and specify font size
        xlab("Position") + ylab("Frequency") + ggtitle(title) + ylim(0, 2.5) + theme(text = element_text(size=14)) +
    
        #Specify increased distance between axes and axes lables so easier to read
        theme(axis.title.y=element_text(margin=margin(0,20,0,0))) + theme(axis.title.x=element_text(margin=margin(20,0,0,0)))

    #Save plot as PDF file
    ggsave(filename=filename, plot=myplot)
}
