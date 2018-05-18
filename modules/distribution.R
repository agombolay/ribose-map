#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP coverage at each chromosome position
#2. Saves plots as png files to appropriate directory

#############################################################################################################################
#Load config
source(commandArgs(TRUE)[1])

#Load libraries
library(ggplot2); library(tools)

#############################################################################################################################
#Specify output directory and file
output <- file.path(repository, "results", sample, "distribution")
input_files <- list.files(path=output, pattern=".tab", full.names=T, recursive=F)
	
for(file in input_files){
        
	#Check size of file > 0
	if (file.info(file)$size > 0){
	
        	#Specify dataset
		data=read.table(file, sep="\t", header=FALSE)
	
		#Re-order levels of strand
		data$V5_new = factor(data$V5, levels=c('+','-'))
	
		#Strand labels for facet wrap
		labels <- c('+' = 'Forward Strand', '-' = 'Reverse Strand')
	
#############################################################################################################################
		#Create plot
		myplot <- ggplot(data, aes(V3, V4, colour = V5_new)) + geom_bar(stat = "identity") +
		
			  #Plot forward and reverse strands, no legend
			  facet_wrap(~V5_new, ncol = 1, labeller = labeller(V5_new = labels)) + 
		
			  #Decrease space between data and axes
			  scale_y_continuous(expand = c(0.015, 0)) + scale_x_continuous(expand = c(0.015, 0)) +
		
			  #Plot scatterplot and set font size
			  theme(text = element_text(size = 20)) + theme(legend.position = "none") + ylim(1, max(data$V4)) +
	
			  #Specify colors for scatterplot and titles for axes
			  xlab("Chromosome Coordinate") + ylab("rNMP Coverage (RPM)") + scale_colour_manual(values = c("blue", "green3"))
	
			  #Remove and replace default background theme of plot
			  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line=element_line(colour = "black"))
		
#############################################################################################################################
		#Save plot as PNG file
		ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=15)

}
}
