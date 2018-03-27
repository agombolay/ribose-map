#!/usr/bin/env Rscript

#Copyright 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP nt frequencies for mito and nucleus
#2. Saves plots as png files to appropriate directory

#############################################################################################################################
#Load config
source(commandArgs(TRUE)[1])

#Load libraries
library(ggplot2); library(tools)
#############################################################################################################################
#Specify output directory and file
output <- file.path(repository, "results", sample, "frequencies")
input_files <- list.files(path=output, pattern=".tab", full.names=T, recursive=F)

for(file in input_files){
	
	#Check size of file > 0
	if (file.info(file)$size > 0){
	
		#Plot regular and zoomed datasets
		for(i in c("normal", "zoomed")) {

#############################################################################################################################
			#Specify datasets to be used for each round of loop
			if (i=="normal") {data=read.table(file, sep="\t", header=TRUE)}
			if (i=="zoomed") {data=read.table(file, sep="\t", header=TRUE)[86:116,]}
    
			#Define variables to store nucleotide positions and frequency values
			position <- data$X; A <- data$A; C <- data$C; G <- data$G; T <- data$U.T

#############################################################################################################################
			#Create plot
			myplot <- ggplot(data, aes(x=position)) +
			
			#Add axes titles and plot title
			xlab("Chromosome Position") + ylab("rNMP Frequency") +
		
			#Specify font size and format legend
			theme(text=element_text(size=20)) + theme(legend.key = element_blank()) +
		
			#Specify color values and remove legend title
                	scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="") +
	
			#Plot data as scatterplot with connecting lines
                	geom_line(aes(y=A, colour="A")) + geom_point(aes(y=A, colour="A")) + geom_line(aes(y=C, colour="C")) + 
			geom_point(aes(y=C, colour="C")) + geom_line(aes(y=G, colour="G")) + geom_point(aes(y=G, colour="G")) +
                	geom_line(aes(y=T, colour="U/T")) + geom_point(aes(y=T, colour="U/T")) +
	
			#Remove and replace default background theme of plot
			theme(panel.grid=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black")) +
			
			#Specify size of legend symbols and remove line through them
			guides(colour=guide_legend(override.aes=list(size=5, linetype=0)))

#############################################################################################################################
			#Save plot as PNG file
			ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".png", sep="")), plot=myplot)
			
}
}
}
