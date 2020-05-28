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
input_files <- list.files(path = output, pattern = "*normalized.tab", full.names = TRUE, recursive = FALSE)	

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


####################################################################################################################################################################
	for(file in input_files){
	
		#Regular and zoomed datasets
		for(i in c("normal", "zoomed")) {

####################################################################################################################################################################
			#Specify datasets to be used for each round of loop
			#if (i == "normal") {data = read.table(file, sep = "\t", header = TRUE)}
			if (i == "zoomed") {data = read.table(file, sep = "\t", header = TRUE)[96:106,]}
    
			#Define variables to store nucleotide positions and frequencies
			position <- data$X; A <- data$A; C <- data$C; G <- data$G; T <- data$U.T

####################################################################################################################################################################
			if (ymin > 0) {

				log <- ggplot(data, aes(x = position)) + theme_minimal() + xlab("Position relative to rNMP") + ylab("Normalized Frequency (log2)") +

				       #Specify color and no legend title
				       scale_colour_manual(values = c("A" = "red2", "C" = "blue4", "G" = "darkorange2", "U/T" = "green4"), name = "") +

				       #Plot data as scatterplot with lines
				       geom_point(aes(y = A, colour = "A"), size = 4, shape = 15) + geom_point(aes(y = C, colour = "C"), size = 4, shape = 16) +
				       geom_point(aes(y = G, colour = "G"), size = 4, shape = 17) + geom_point(aes(y = T, colour = "U/T"), size = 4, shape = 18) +
				       geom_line(aes(y = A, colour = "A"), size = .75) + geom_line(aes(y = C, colour = "C"), size = .75) +
				       geom_line(aes(y = G, colour = "G"), size = .75) + geom_line(aes(y = T, colour = "U/T"), size = .75) +		   

				       #Format legend symbols and specify y-axis limits
				       guides(colour = guide_legend(override.aes = list(size = 5, linetype = 0, shape = c(15, 16, 17, 18)))) + scale_y_continuous(trans = 'log2', limits = c(ymin, ymax), labels = scales::number_format(accuracy = 0.01)) +
				       scale_x_continuous(limits=c(-5,5)) +
				
				       #Add axis lines and ticks and increase font size
				       theme(
					     axis.title = element_text(color = "black", size = 25), axis.line = element_line(size = 1), axis.text = element_text(color = "black", size = 25),
				  	     axis.ticks = element_line(colour = "black", size = 1), axis.ticks.length = unit(.4, "cm"),legend.text = element_text(color = "black", size = 20),
				  	     axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
				  	     plot.margin = unit(c(.5, .5, .5, .5), "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()
				       )


			        ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".png", sep = "")), plot = log)
				ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".pdf", sep = "")), plot = log)
			
####################################################################################################################################################################
			} else {
				linear <- ggplot(data, aes(x = position)) + theme_minimal() + xlab("Position relative to rNMP") + ylab("Normalized Frequency") +

				       #Specify color and no legend title
				       scale_colour_manual(values = c("A" = "red2", "C" = "blue4", "G" = "darkorange2", "U/T" = "green4"), name = "") +

				       #Plot data as scatterplot with lines
				       geom_point(aes(y = A, colour = "A"), size = 4, shape = 15) + geom_point(aes(y = C, colour = "C"), size = 4, shape = 16) +
				       geom_point(aes(y = G, colour = "G"), size = 4, shape = 17) + geom_point(aes(y = T, colour = "U/T"), size = 4, shape = 18) +
				       geom_line(aes(y = A, colour = "A"), size = .75) + geom_line(aes(y = C, colour = "C"), size = .75) +
				       geom_line(aes(y = G, colour = "G"), size = .75) + geom_line(aes(y = T, colour = "U/T"), size = .75) +		   

				       #Format legend symbols and specify y-axis limits
				       guides(colour = guide_legend(override.aes = list(size = 5, linetype = 0, shape = c(15, 16, 17, 18)))) + scale_y_continuous(limits = c(0, ymax)) +

				       #Add axis lines and ticks and increase font size
				       theme(
					     axis.title = element_text(color = "black", size = 25), axis.line = element_line(size = 1), axis.text = element_text(color = "black", size = 25),
				  	     axis.ticks = element_line(colour = "black", size = 1), axis.ticks.length = unit(.4, "cm"), legend.text = element_text(color = "black", size = 20),
				  	     axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
				  	     plot.margin = unit(c(.5, .5, .5, .5), "cm")
				       )

			        ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".png", sep = "")), plot = linear)
				ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".pdf", sep = "")), plot = linear)
}
}
}
message("Status: Sequence Module plotting for ", sample, " is complete")
warnings()
