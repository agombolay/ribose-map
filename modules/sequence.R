#!/usr/bin/env Rscript

#Copyright 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP nt frequencies for mito and nucleus
#2. Saves plots as png files to appropriate directory

####################################################################################################################################################################
#Load config
source(commandArgs(TRUE)[1])

#Load libraries
library(ggplot2); library(tools)

####################################################################################################################################################################
#Input/Output
output <- file.path(repository, "results", sample, "sequence")
input_files <- list.files(path = output, pattern = ".tab", full.names = TRUE, recursive = FALSE)

####################################################################################################################################################################
#Find maximum y-axis value
maximum <- c()

for(file in input_files){
	if (file.info(file)$size > 0){
		data = read.table(file, sep = "\t", header = TRUE)
		A <- data$A; C <- data$C; G <- data$G; T <- data$U.T
		maximum <- c(maximum, max(A, C, G, T))
}
}
ylimit <- max(maximum)

####################################################################################################################################################################
for(file in input_files){
	
	#Check size of file > 0
	if (file.info(file)$size > 0){
	
		#Regular and zoomed datasets
		for(i in c("normal", "zoomed")) {

####################################################################################################################################################################
			#Specify datasets to be used for each round of loop
			if (i == "normal") {data = read.table(file, sep = "\t", header = TRUE)}
			if (i == "zoomed") {data = read.table(file, sep = "\t", header = TRUE)[86:116,]}
    
			#Define variables to store nucleotide positions and frequencies
			position <- data$X; A <- data$A; C <- data$C; G <- data$G; T <- data$U.T

####################################################################################################################################################################
			combined <- ggplot(data, aes(x = position)) + xlab("Chromosome Position") + ylab("Normalized Nucleotide Frequency") +
			
					   #Specify color and remove legend title
                			   scale_colour_manual(values = c("A" = "#CC79A7", "C" = "#56B4E9", "G" = "#E69F00", "U/T" = "#009E73"), name = "") +

					   #Plot as scatterplot with connecting lines
					   geom_point(aes(y = A, colour = "A")) + geom_point(aes(y = C, colour = "C")) + geom_point(aes(y = G, colour = "G")) +
					   geom_point(aes(y = T, colour = "U/T")) + geom_line(aes(y = A, colour = "A")) + geom_line(aes(y = C, colour = "C")) +
					   geom_line(aes(y = G, colour = "G")) + geom_line(aes(y = T, colour = "U/T")) +

					   #Format legend symbols (no line through symbol and increase symbol size)
					   guides(colour = guide_legend(override.aes = list(size = 5, linetype = 0))) + theme_minimal() +
			
					   #Simplify background theme, increase font size, and specify y-axis limits
					   theme(text = element_text(size = 20)) + scale_y_continuous(limits = c(0, ylimit)) +
					   
					   theme(axis.line = element_line(color="black", size = 1))
			
					   #Replace default background with black lines for the axes
					   #annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) + annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)
		
####################################################################################################################################################################
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", "combined", ".", i, ".png", sep = "")), plot = combined)
			
}
}
}
