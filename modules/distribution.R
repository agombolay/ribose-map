#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP coverage at each chromosome position
#2. Saves plots as png files to appropriate directory

####################################################################################################################################################################
#Load config
source(commandArgs(TRUE)[1])

#Load libraries
library(ggplot2); library(tools)

####################################################################################################################################################################
#Input/Output
output <- file.path(repository, "results", sample, paste("distribution-", quality, sep=""))
input_files <- list.files(path = output, pattern = ".tab", full.names = TRUE, recursive = FALSE)

####################################################################################################################################################################
#Find maximum y-axis value
maximum <- c()

for(file in input_files){
	if (file.info(file)$size > 0){
		data = read.table(file, sep = "\t", header = FALSE)
		maximum <- c(maximum, max(data$V5))
}
}
ylimit <- max(maximum)

####################################################################################################################################################################
for(file in input_files){
        
	#Check size of file > 0
	if (file.info(file)$size > 0){
	
        	#Specify dataset
		data = read.table(file, sep = "\t", header = FALSE)
		
		#Re-order levels of DNA strand
		data$V4_new = factor(data$V4, levels = c('+', '-'))
		
		#Specify DNA strand labels for plot
		labels <- c('+' = 'Forward Strand', '-' = 'Reverse Strand')
		
####################################################################################################################################################################
		if (scale == 'yes'){
			
			scaled <- ggplot(data, aes(V3, V5, colour = V4_new)) + xlab("Chromosome Coordinate") + ylab("Per Nucleotide rNMP Coverage (%)") +

				 #Plot forward and reverse strands on same plot
				 facet_wrap(~V4_new, ncol = 1, labeller = labeller(V4_new = labels)) + scale_colour_manual(values = c("#0072B2", "#009E73")) +

				 #Plot barcharts, specify font size and no legend, and remove background
				 geom_bar(stat = "identity") + theme(text = element_text(size = 20), legend.position = "none", panel.background = element_blank()) +		 

				 #Add lines to represent axes (must be added for plots that use facet wrap)
				 annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) + annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +

				 #Decrease space between barcharts and x/y-axes
				 scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0), limits = c(0, ylimit))

			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep = "")), plot = scaled, width = 15)
		
####################################################################################################################################################################
		} else if (scale == 'no'){
	
			unscaled <- ggplot(data, aes(V3, V5, colour = V4_new)) + xlab("Chromosome Coordinate") + ylab("Per Nucleotide rNMP Coverage (%)") +

				 #Plot forward and reverse strands on same plot
				 facet_wrap(~V4_new, ncol = 1, labeller = labeller(V4_new = labels)) + scale_colour_manual(values = c("#0072B2", "#009E73")) +

				 #Plot barcharts, specify font size and no legend, and remove background
				 geom_bar(stat = "identity") + theme(text = element_text(size = 20), legend.position = "none", panel.background = element_blank()) +		 

				 #Add lines to represent axes (must be added for plots that use facet wrap)
				 annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) + annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +

				 #Decrease space between barcharts and x/y-axes
				 scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0))

			ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep = "")), plot = unscaled, width = 15)

}
}
}

message("Status: Distribution Module plotting for ", sample, " is complete")
