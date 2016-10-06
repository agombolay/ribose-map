#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
title <- args[2]

#Define data file to read
data <- read.table(file, sep="\t", header=TRUE)

#Define data values
position <- data$X
frequencyA <- data$A
frequencyC <- data$C
frequencyG <- data$G
frequencyT <- data$U.T

#Load ggplot2 package
library(ggplot2)

#Plot data values using ggplot2
ggplot(data=data, aes(x=position)) +
    geom_line(aes(y = frequencyA, colour = "A")) + geom_point(aes(y = frequencyA, colour = "A")) +
    geom_line(aes(y = frequencyC, colour = "C")) + geom_point(aes(y = frequencyC, colour = "C")) +
    geom_line(aes(y = frequencyG, colour = "G")) + geom_point(aes(y = frequencyG, colour = "G")) +
    geom_line(aes(y = frequencyT, colour = "U/T")) + geom_point(aes(y = frequencyT, colour = "U/T")) +
    xlab("Position Upstream (-)/Downstream (+) from rNMP") + ylab("Normalized Frequency") + ggtitle(title) +
    scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="Nucleotide") +  theme_bw() +
    theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + theme(panel.border= element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) + theme(legend.key = element_rect(colour = NA))
