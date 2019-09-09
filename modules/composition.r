#!/usr/bin/env Rscript

#Load config
source(commandArgs(TRUE)[1])

#Load libraries
library(ggplot2); library(tools)

#Input/Output
output <- file.path(repository, "results", sample, paste("composition", quality, sep = ""))
input_files <- list.files(path = output, pattern = "*.tab", full.names = TRUE, recursive = FALSE)

for(file in input_files){

	data = read.table(file, sep = "\t", header = FALSE)

	barplot <- ggplot(data, aes(x = data$V1, y = data$V2, fill = data$V1)) + geom_bar(stat = "identity", width = 0.75) + xlab("") + ylab("%") +

	           scale_fill_manual(values = c("red2", "blue4", "darkorange2", "green4"), name = "") +

	           scale_x_discrete(expand = c(0.2, 0)) + scale_y_continuous(expand = c(0.01, 0)) +

	           theme(
    	             axis.title = element_text(color = "black", size = 25), axis.line = element_line(size = 1), axis.text = element_text(color = "black", size = 25),
    	             axis.ticks = element_line(colour = "black", size = 1), axis.ticks.length = unit(.4, "cm"), legend.text = element_text(color = "black", size = 20),
    	             axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    	             plot.margin = unit(c(.5, .5, .5, .5), "cm"), panel.background = element_blank()
		       )
    ggsave(filename = file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".png", sep = "")), plot = barplot)

}
