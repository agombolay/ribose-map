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
output <- file.path(repository, "results", sample, "sequence")
input_files <- list.files(path = output, pattern = ".tab", full.names = T, recursive = F)

for(file in input_files){
	
	#Check size of file > 0
	if (file.info(file)$size > 0){
	
		#Plot regular and zoomed datasets
		for(i in c("normal", "zoomed")) {

#############################################################################################################################
			#Specify datasets to be used for each round of loop
			if (i == "normal") {data = read.table(file, sep = "\t", header = TRUE)}
			if (i == "zoomed") {data = read.table(file, sep = "\t", header = TRUE)[86:116,]}
    
			#Define variables to store nucleotide positions and frequency values
			position <- data$X; A <- data$A; C <- data$C; G <- data$G; T <- data$U.T

#############################################################################################################################
			#Create plot
			combined <- ggplot(data, aes(x = position)) + xlab("Chromosome Position") + ylab("Nucleotide Frequency") +
		
				           #Specify color and remove legend title
                		    	   scale_colour_manual(values = c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name = "") +
	
				    	   #Plot as scatterplot with connecting lines
                			   geom_line(aes(y = A, colour = "A")) + geom_point(aes(y = A, colour = "A")) + geom_line(aes(y = C, colour = "C")) + 
					   geom_point(aes(y = C, colour = "C")) + geom_line(aes(y = G, colour = "G")) + geom_point(aes(y = G, colour = "G")) +
                			   geom_line(aes(y = T, colour = "U/T")) + geom_point(aes(y = T, colour = "U/T")) +
			
					   #Format legend symbols and specify font size
					   guides(colour = guide_legend(override.aes = list(size = 5, linetype = 0))) + theme(text = element_text(size = 20)) +
	
					   #Simplify default ggplot2 background formatting
				  	   theme(legend.key = element_blank()) + theme(panel.background = element_blank(), axis.line=element_line(colour = "black")) +
					   ylim(0, max(A,C,G,T))
			
#############################################################################################################################
			nucleotideA <- ggplot(data, aes(x = position)) + xlab("Chromosome Position") + ylab("Nucleotide Frequency") + theme(legend.position = "none") +
		
					      theme(text = element_text(size = 20)) + scale_colour_manual(values = c("#CC79A7"), name="") +
                			      geom_line(aes(y = A, colour = "A")) + geom_point(aes(y = A, colour = "A")) +
					      theme(panel.background=element_blank(), axis.line = element_line(colour = "black")) + ylim(0, max(A,C,G,T))

			nucleotideC <- ggplot(data, aes(x = position)) + xlab("Chromosome Position") + ylab("Nucleotide Frequency") + theme(legend.position = "none") +
		
					      theme(text = element_text(size = 20)) + scale_colour_manual(values = c("#56B4E9"), name = "") +
                			      geom_line(aes(y = C, colour = "C")) + geom_point(aes(y = C, colour = "C")) +
					      theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
			
			nucleotideG <- ggplot(data, aes(x = position)) + xlab("Chromosome Position") + ylab("Nucleotide Frequency") + theme(legend.position = "none") +
		
					      theme(text = element_text(size=20)) + scale_colour_manual(values = c("#E69F00"), name = "") +
                			      geom_line(aes(y = G, colour = "G")) + geom_point(aes(y = G, colour = "G")) +
					      theme(panel.background = element_blank(), axis.line=element_line(colour = "black"))

			nucleotideT <- ggplot(data, aes(x = position)) + xlab("Chromosome Position") + ylab("Nucleotide Frequency") + theme(legend.position = "none") +
		
                			      geom_line(aes(y = T, colour = "T")) + geom_point(aes(y = T, colour = "T")) + scale_colour_manual(values = c("#009E73"), name = "") +
					      theme(text = element_text(size = 20)) + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

#############################################################################################################################
			#Save plot as PNG file
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", "combined", ".", i, ".png", sep = "")), plot = combined)
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", "nucleotideA", ".", i, ".png", sep = "")), plot = nucleotideA)
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", "nucleotideC", ".", i, ".png", sep = "")), plot = nucleotideC)
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", "nucleotideG", ".", i, ".png", sep = "")), plot = nucleotideG)
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", "nucleotideT", ".", i, ".png", sep = "")), plot = nucleotideT)
}
}
}
