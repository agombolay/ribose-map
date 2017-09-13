#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots the number of rNMPs at each genomic position.

#Load libraries
library(optparse)
library(ggplot2)
library(tools)

#Command line options
option_list <- list(
make_option(c("-s", "--sample"), help="Sample name(s) (e.g., FS1, FS2, FS3)"),
make_option(c("-t", "--title"), help="Title will be same for all plots; blank space = no titles"),
make_option(c("-r", "--reference"), help="Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)"),
make_option(c("-d", "--directory"), help="Directory (e.g., /projects/home/agombolay3/data/repository/Ribose-Map)")
)

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

for(i in opt$sample) {
        
        #Specify output directory and file
        output <- file.path(opt$directory, "Results", opt$reference, opt$sample, "Hotspots")
	files <- list.files(path=output, pattern=".bed", full.names=TRUE, recursive=FALSE)
        
	for(file in files){
		
		#Plot only if files exist
        	if (file.exists(file)) {
            
        		#Specify dataset
            		data=read.table(file, sep="\t", header=FALSE)

#############################################################################################################################
        		#Plot coverage
        		myplot <- ggplot(data, aes(x=1:length(data[,3]), y=data[,3])) +

			#Replace default theme
                	theme(panel.grid=element_blank(),
                      	      panel.background=element_blank(),
                      	      axis.line=element_line(colour="black")) +
			
        		#Add axes titles and plot title
        		xlab("Position in genome") + ylab("Frequency of rNMPs") + ggtitle("") +

			#Plot data as scatterchart with connecting lines
        		geom_point(shape=1, colour="blue4") + geom_line(aes(y=data[,3]), colour="blue4") +
				
        		#Decrease space between scatterplot and x-axis/y-axis
        		scale_y_continuous(expand=c(0.015,0)) + scale_x_continuous(expand=c(0.015,0)) +

        		#Specify font size and center title (if any) of plot on page
        		theme(text=element_text(size=20)) + theme(plot.title=element_text(hjust=0.5))

#############################################################################################################################
#Save plot as PNG file
ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=20)

}
}
}
