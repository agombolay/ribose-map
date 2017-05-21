#!/usr/bin/env Rscript

#Copyright 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots rNMP and flanking nucleotide frequencies (100 bp upstream and downstream)
#rNMP frequencies = position 0, upstream frequencies = -100 -> -1, and downstream = +1 -> +100

#Load libraries
library(optparse)
library(ggplot2)

#Command line options
option_list <- list(
    make_option(c("-i", "--input")),
    make_option(c("-t", "--title")),
    make_option(c("-s", "--sample"))
)

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

#Plot regular dataset
if (file.exists(opt$input)) {
    #Define data file (tab delimited) to read
    dataset <- read.table(opt$input, sep="\t", header=TRUE)
    
    #Define position values
    position <- dataset$X
    
    #Define frequency values
    frequencyA <- dataset$A; frequencyC <- dataset$C
    frequencyG <- dataset$G; frequencyT <- dataset$U.T

    #Plot frequencies
    myplot <- ggplot(data=dataset, aes(x=position)) +
        
        #Plot data as scatterplot   
        geom_line(aes(y = frequencyA, colour = "A")) + geom_point(aes(y = frequencyA, colour = "A")) +
        geom_line(aes(y = frequencyC, colour = "C")) + geom_point(aes(y = frequencyC, colour = "C")) +
        geom_line(aes(y = frequencyG, colour = "G")) + geom_point(aes(y = frequencyG, colour = "G")) +
        geom_line(aes(y = frequencyT, colour = "U/T")) + geom_point(aes(y = frequencyT, colour = "U/T")) +
    
        #Add axes titles
        xlab("Position") + ylab("Frequency") +

        #Remove and replace default background plot theme
        theme_bw() + theme(panel.border=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) +
    
        #Specify color values for each of the four nucleotides
        scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="Nucleotide") +
    
        #Label axes, add title, define y-axis limits, and specify font size
        theme(text = element_text(size=14)) + ggtitle(opt$title) + theme(plot.title = element_text(hjust = 0.5)) +

        #Specify increased distance between axes and axes lables so easier to read
        theme(axis.title.y=element_text(margin=margin(0,20,0,0))) + theme(axis.title.x=element_text(margin=margin(20,0,0,0)))

    #Save plot as PDF file
    ggsave(filename=paste(opt$sample, "-Regular", ".png", sep=""), plot=myplot)
}

#Plot zoomed dataset
if (file.exists(opt$input)) {
    #Define data file (tab delimited) to read
    zoomed <- dataset[86:116,]
    
    #Define position values
    position <- zoomed$X
    
    #Define frequency values
    frequencyA <- zoomed$A; frequencyC <- zoomed$C
    frequencyG <- zoomed$G; frequencyT <- zoomed$U.T
    
    #Plot frequencies
    myplot <- ggplot(data=zoomed, aes(x=position)) +
    
    #Plot data as scatterplot
    geom_line(aes(y = frequencyA, colour = "A")) + geom_point(aes(y = frequencyA, colour = "A")) +
    geom_line(aes(y = frequencyC, colour = "C")) + geom_point(aes(y = frequencyC, colour = "C")) +
    geom_line(aes(y = frequencyG, colour = "G")) + geom_point(aes(y = frequencyG, colour = "G")) +
    geom_line(aes(y = frequencyT, colour = "U/T")) + geom_point(aes(y = frequencyT, colour = "U/T")) +
    
    #Add axes titles
    xlab("Position") + ylab("Frequency") +
    
    #Remove and replace default background plot theme
    theme_bw() + theme(panel.border=element_blank(), panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) +
    
    #Specify color values for each of the four nucleotides
    scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="Nucleotide") +
    
    #Label axes, add title, define y-axis limits, and specify font size
    theme(text = element_text(size=14)) + ggtitle(opt$title) + theme(plot.title = element_text(hjust = 0.5)) +
    
    #Specify increased distance between axes and axes lables so easier to read
    theme(axis.title.y=element_text(margin=margin(0,20,0,0))) + theme(axis.title.x=element_text(margin=margin(20,0,0,0)))
    
    #Save plot as PDF file
    ggsave(filename=paste(opt$sample, "-Zoomed", ".png", sep=""), plot=myplot)
}
