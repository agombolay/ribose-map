#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP coverage at each chromosome position
#2. Saves plots as png files to appropriate directory

#############################################################################################################################
#Load library
library(ggplot2)

#Load config file
source(commandArgs(TRUE)[1])

#############################################################################################################################
#Specify output directory and file
output <- file.path(directory, "results", sample, "distribution")
input_files <- list.files(path=output, pattern=".bed", full.names=T, recursive=F)
        
for(file in input_files){
            
        #Specify dataset
	data=read.table(file, sep="\t", header=FALSE)
			
	labels <- c('+' = 'Forward', '-' = 'Reverse')
	data$V5_new = factor(data$V5, levels=c('+','-'))
	
#############################################################################################################################
	#Create plot
	myplot <- ggplot(data,aes(V3, V4, colour=V5_new)) + geom_point() + ylim(1,max(data$V4)) +
	
	#Plot forward and reverse strands, no legend
	facet_wrap(~V5_new,ncol=1,labeller=labeller(V5_new = labels)) + theme(legend.position="none") +
	
	#Specify colors for scatterplot and titles for axes
	scale_colour_manual(values=c("blue", "green3")) + xlab("Chromosome Position") + ylab("rNMP Frequency") +
	
	#Remove and replace default background theme of plot
	theme(panel.grid=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"))

#############################################################################################################################
	#Save plot as PNG file
	ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=20)

}
