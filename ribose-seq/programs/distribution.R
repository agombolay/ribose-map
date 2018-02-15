#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP coverage at each chromosome position
#2. Saves plots as png files to appropriate directory

#############################################################################################################################
#Load config file
path <- commandArgs(TRUE)[1]
source(file.path(path,"config.R"))

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
	#myplot <- ggplot(data, aes(x=data[,3], y=values)) +
		
	myplot <- ggplot(data, aes(x=data[,3], y=values)) +
        geom_segment(aes(x = 0, y = values, xend = data[,3], yend = values), color = "grey50") +
        geom_point()
	
	#Remove and replace default background theme of plot
	#theme(panel.grid=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black")) +
	
	#Create barchart and add titles for axes
        #geom_bar(stat="identity", fill="black", color="black", width = .0001) + xlab("Chromosome Position") + ylab("rNMP Frequency")+ 
	
	#Decrease spacing between plot and axes and increase font
	#theme(text=element_text(size=15))

#############################################################################################################################
	#Save plot as PNG file
	ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=20)

}
