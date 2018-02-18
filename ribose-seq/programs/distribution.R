#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP coverage at each chromosome position
#2. Saves plots as png files to appropriate directory

#############################################################################################################################
#Load config file
config <- commandArgs(TRUE)[1]
#source(file.path(path,"config.R"))
source(config)

#Load R libraries
library(tools); library(ggplot2); library(optparse)

#############################################################################################################################
#Specify output directory and file
output <- file.path(directory, "results", sample, "distribution")
input_files <- list.files(path=output, pattern=".bed", full.names=T, recursive=F)
        
for(file in input_files){
            
        #Specify dataset
	data=read.table(file, sep="\t", header=FALSE)
			
	#Transform values on negative strand
	forward <- data[ which(data$V5=='+'), ]; reverse <- data[ which(data$V5=='-'), ]

#############################################################################################################################
	#Create plot
	myplot <- ggplot() + geom_point(aes(x=forward$V3, y=forward$V4, colour="forward")) +
	
	geom_point(aes(x=reverse$V3, y=reverse$V4*-1, colour="reverse")) + theme(text=element_text(size=20))
	
	scale_colour_manual(values=c("#CC79A7", "#56B4E9"), name="") + theme(legend.key = element_blank()) +
	
        xlab("Chromosome Position") + ylab("rNMP Frequency")
	
	#Remove and replace default background theme of plot
	theme(panel.grid=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"))

#############################################################################################################################
	#Save plot as PNG file
	ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=20)

}
