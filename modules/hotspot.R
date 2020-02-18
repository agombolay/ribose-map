#!/usr/bin/env Rscript

#Load the required packages
library(tools); require(ggplot2); require(ggseqlogo)

#Input/Output
input <- list.files(path = "/nv/hp16/agombolay3/data/agombolay3_copy/algae/consensus/percentile_90/alternative_plotting", pattern = "*.txt", full.names = TRUE, recursive = TRUE)

for(file in input){

	if (file.size(file) > 0) {

	#Some sample data
	data <- readLines(file)

	#Select sequences
	rows <- grep("\\s\\.\\s", data, value = TRUE)
	sequences <- substr(rows, 54, 60)

	#Specify color scheme
	colors = make_col_scheme(chars = c('A', 'C', 'G', 'T'), cols = c('red2', 'blue4', 'darkorange2', 'green4'))

	#Plot customized logo
	plot <- ggseqlogo(sequences, seq_type = 'dna', col_scheme = colors, method = 'bits') + theme_classic() + xlab("") + ylab("") +

	scale_y_continuous(limits = c(0, 0.5), expand = c(0.02, 0)) + scale_x_continuous(breaks = seq(1,7,1), labels = seq(-3,3,1), expand = c(0.02, 0)) +

	theme(
		  axis.line = element_line(color = "black", size = 1), axis.text = element_text(color = "black", size = 45), axis.ticks = element_line(color
		  = "black", size = 1), axis.ticks.length = unit(0.5, "cm"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
	)

	ggsave(filename = file.path(paste(dirname(file),"-","ggseqlogo",".png",sep = "")), plot = plot)
}
}
