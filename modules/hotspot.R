#!/usr/bin/env Rscript

#Load config
source(commandArgs(TRUE)[1])

#Load the required packages
library(tools); require(ggplot2); require(ggseqlogo)

#Input/Output
output <- file.path(repository, "results", sample, paste("hotspot", quality, sep = ""))
input_files <- list.files(path = output, pattern = "*.txt", full.names = TRUE, recursive = TRUE)

for(file in input_files){

	if (file.size(file) > 0) {

	#Some sample data
	data <- readLines(file)

	#Select sequences
	rows <- grep("\\s\\.\\s", data, value = TRUE)
	sequences <- substr(rows, 54, 60)

	#Specify color scheme
	colors = make_col_scheme(chars = c('A', 'C', 'G', 'T'), cols = c('red2', 'blue4', 'darkorange2', 'green4'))

	#Plot customized logo
	plot <- ggseqlogo(sequences, seq_type = 'dna', col_scheme = colors, method = 'bits') + theme_classic() + xlab("Position relative to rNMP") + ylab("Bits") +

	scale_y_continuous(limits = c(0, 2), expand = c(0.02, 0)) + scale_x_continuous(breaks = seq(1,7,1), labels = seq(-3,3,1), expand = c(0.02, 0)) +

	theme(
	      axis.title = element_text(color = "black", size = 25), axis.line = element_line(color = "black", size = 1), axis.text = element_text(color = "black", size = 25),
	      axis.ticks = element_line(color = "black", size = 1), axis.ticks.length = unit(0.5, "cm"), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
	      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
	)

	ggsave(filename = file.path(paste(dirname(file),"-","ggseqlogo",".png",sep = "")), plot = plot)
}
}
