#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates the Chi Square statistic and its associated p-value.
#If the observed or expected counts are < 5, then categories will be collapsed.

#Load libraries
library(optparse)
library(ggplot2)

#Command line options
option_list <- list(
make_option(c("-s", "--sample"), help="Sample name(s) (e.g., FS1, FS2, FS3)"),
make_option(c("-t", "--title"), help="Title will be same for all plots; blank space = no titles"),
make_option(c("-r", "--reference"), help="Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)"),
make_option(c("-d", "--directory"), help="Local user directory (e.g., /projects/home/agombolay3/data/repository)")
)

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

for(i in opt$sample) {
    for(j in c("nucleus")) {
        
        #Specify output directory and file
        output <- file.path(opt$directory, "Ribose-Map", "Results", opt$reference, opt$sample, "Hotspots")
        file <- file.path(output, paste(opt$sample, "-", "Coverage", ".", opt$reference, ".", j, ".txt", sep=""))

        #Plot only if files exist
        if (file.exists(file)) {
            
            #Specify dataset
            data=read.table(file, sep="\t", header=TRUE)

#############################################################################################################################
        #Plot coverage
        myplot <- ggplot(data, aes(x=1:length(data[,3]), y=data[,3])) +

        #Add axes titles and plot title
        xlab("Position in genome") + ylab("Frequency of rNMPs") + ggtitle("") +

        #Replace default background plot theme
        theme_bw() + theme(panel.border=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) +

        #Plot data as scatterchart with connecting lines
        geom_point(shape=1, colour="blue4") + geom_line(aes(y=data[,3]), colour="blue4") +
				
        #Decrease space between scatterplot and x-axis/y-axis
        scale_y_continuous(expand=c(0.015,0)) + scale_x_continuous(expand=c(0.015,0)) +

        #Specify font size and center title (if any) of plot on page
        theme(text=element_text(size=14)) + theme(plot.title=element_text(hjust=0.5))

#############################################################################################################################
#Save plot as PNG file
ggsave(filename=file.path(output, paste(opt$sample, "-", "Hotspots", ".", j, ".png", sep="")), plot=myplot, height=10, width=12)

message("Plotting of", " ", i , " ", "is complete")
}
}
}
