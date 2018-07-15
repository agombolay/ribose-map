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
output <- file.path(repository, "results", sample, "sequence", "-", quality)

if (plots == 'all'){
	input_files <- list.files(path = output, pattern = "*.tab", full.names = TRUE, recursive = FALSE)
} else if (plots == 'combined'){
	input_files <- list.files(path = output, pattern = "*Combined.*.tab", full.names = TRUE, recursive = FALSE)
}

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
print(ylimit)
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
			myplot <- ggplot(data, aes(x = position)) + theme_minimal() + xlab("Chromosome Position") + ylab("Normalized Frequency") +

					   #Specify color and no legend title
                			   scale_colour_manual(values = c("A" = "#CC79A7", "C" = "#56B4E9", "G" = "#E69F00", "U/T" = "#009E73"), name = "") +

					   #Plot data as scatterplot with lines
					   geom_point(aes(y = A, colour = "A")) + geom_point(aes(y = C, colour = "C")) + geom_point(aes(y = G, colour = "G")) +
					   geom_point(aes(y = T, colour = "U/T")) + geom_line(aes(y = A, colour = "A")) + geom_line(aes(y = C, colour = "C")) +
					   geom_line(aes(y = G, colour = "G")) + geom_line(aes(y = T, colour = "U/T")) +		   

					   #Format legend symbols and specify y-axis limits
					   guides(colour = guide_legend(override.aes = list(size = 5, linetype = 0))) + scale_y_continuous(limits = c(0, ylimit)) +

					   #Add axis lines and ticks and increase font size
					   theme(axis.line = element_line(size = .5), text = element_text(size = 20), axis.ticks = element_line(colour = "black"))

####################################################################################################################################################################
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".png", sep = "")), plot = myplot)
			
}
}
}

message("Status: Sequence Module plotting for ", sample, " is complete")
