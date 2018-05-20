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
input_files <- list.files(path = output, pattern = ".tab", full.names = T, recursive = F)

#Find maximum y-axis value
maximum <- c()

for(file in input_files){

	if (file.info(file)$size > 0){
		
		data = read.table(file, sep = "\t", header = F)
		maximum <- c(maximum, max(data$V4))
}
}

ylimit <- max(maximum)

for(file in input_files){
        
	#Check size of file > 0
	if (file.info(file)$size > 0){
	
        	#Specify dataset
		data = read.table(file, sep = "\t", header = F)
		
		#Re-order levels of DNA strand
		data$V5_new = factor(data$V5, levels = c('+','-'))
		
		#Specify DNA strand labels for plot
		labels <- c('+' = 'Forward Strand', '-' = 'Reverse Strand')

#############################################################################################################################
		#Create plot
		myplot <- ggplot(data, aes(V3, V4, colour = V5_new)) + geom_bar(stat = "identity") +
		
				 #Plot forward and reverse strands separately
				 facet_wrap(~V5_new, ncol = 1, labeller = labeller(V5_new = labels)) + 
		
				 #Specify font size, no legend, and y-axis limits
				 theme(text = element_text(size = 20)) + theme(legend.position = "none") +
	
				 #Specify titles for axes and colors for barcharts
				 xlab("Chromosome Coordinate") + ylab("rNMP Coverage (%)") + scale_colour_manual(values = c("blue", "green3")) +
	
				 #Remove and replace default background theme of plot
				 theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line=element_line(colour = "black")) +
		
				 #Specify limits and break points for axis
				 scale_y_continuous(limits = c(0,ylimit))
#############################################################################################################################
		#Save plot as PNG file
		ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep = "")), plot = myplot, width = 15)

}
}
