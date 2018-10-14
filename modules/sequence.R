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
output <- file.path(repository, "results", sample, paste("sequence", quality, sep = ""))

if (plots == 'combined'){
	input_files <- list.files(path = output, pattern = "*(Combined)\\.(nucleus|mitochondria)\\.tab$", full.names = TRUE, recursive = FALSE)	
} else if (plots == 'all'){
	input_files <- list.files(path = output, pattern = "*(A|G|C|T|Combined)\\.(nucleus|mitochondria)\\.tab$", full.names = TRUE, recursive = FALSE)
}

####################################################################################################################################################################
#Find maximum y-axis value
maximum <- c()
minimum <- c()

for(file in input_files){
	data = read.table(file, sep = "\t", header = TRUE)
	A <- data$A; C <- data$C; G <- data$G; T <- data$U.T
	maximum <- c(maximum, max(A, C, G, T, na.rm = TRUE))
	minimum <- c(minimum, min(A, C, G, T, na.rm = TRUE))
}
ymax <- max(maximum)
ymin <- min(minimum)

if (ymin > 0) {
####################################################################################################################################################################
	for(file in input_files){
	
		#Regular and zoomed datasets
		for(i in c("normal", "zoomed")) {

####################################################################################################################################################################
			#Specify datasets to be used for each round of loop
			if (i == "normal") {data = read.table(file, sep = "\t", header = TRUE)}
			if (i == "zoomed") {data = read.table(file, sep = "\t", header = TRUE)[86:116,]}
    
			#Define variables to store nucleotide positions and frequencies
			position <- data$X; A <- data$A; C <- data$C; G <- data$G; T <- data$U.T

####################################################################################################################################################################
			myplot <- ggplot(data, aes(x = position)) + theme_minimal() + xlab("Position relative to rNMP") + ylab("Normalized Frequency (log2)") +

				  #Specify color and no legend title
				  scale_colour_manual(values = c("A" = "#CC79A7", "C" = "#56B4E9", "G" = "#E69F00", "U/T" = "#009E73"), name = "") +

				  #Plot data as scatterplot with lines
				  geom_point(aes(y = A, colour = "A")) + geom_point(aes(y = C, colour = "C")) + geom_point(aes(y = G, colour = "G")) +
				  geom_point(aes(y = T, colour = "U/T")) + geom_line(aes(y = A, colour = "A")) + geom_line(aes(y = C, colour = "C")) +
				  geom_line(aes(y = G, colour = "G")) + geom_line(aes(y = T, colour = "U/T")) +		   

				  #Format legend symbols and specify y-axis limits
				  guides(colour = guide_legend(override.aes = list(size = 5, linetype = 0))) + scale_y_continuous(trans = 'log2', breaks = c(ymin, 1, ymax), limits = c(ymin, ymax)) +

				  #Add axis lines and ticks and increase font size
				  theme(
					axis.title = element_text(color = "black", size = 25), axis.line = element_line(size = 1), axis.text = element_text(color = "black", size = 25),
				  	axis.ticks = element_line(colour = "black", size = 1), axis.ticks.length = unit(.4, "cm"),legend.text = element_text(color = "black", size = 20),
				  	axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
				  	plot.margin = unit(c(.5, .5, .5, .5), "cm")
				  )

####################################################################################################################################################################
			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".png", sep = "")), plot = myplot)
			
}
}
}
message("Status: Sequence Module plotting for ", sample, " is complete")
