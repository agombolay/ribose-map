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
	values <- ifelse(data$V5=='-',data$V4*-1,data$V4)

#############################################################################################################################
	#Create plot
	myplot <- ggplot(data, aes(x=data[,3], y=values)) + geom_point() +
	
	#scale_y_discrete(expand=c(0.015,0)) + scale_x_discrete(expand=c(0.015,0)) +
	
	#Remove and replace default background theme of plot
	theme(panel.grid=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black")) +
	
	#Decrease spacing between plot and axes and increase font
        xlab("Chromosome Position") + ylab("rNMP Frequency") + theme(text=element_text(size=15)) +
	
	#Specify color values and remove legend title
	scale_colour_manual(values=c("#CC79A7", "#56B4E9"), name="")

#############################################################################################################################
	#Save plot as PNG file
	ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=20)

}
