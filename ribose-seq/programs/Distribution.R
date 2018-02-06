#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots the number of rNMPs at each chromosome position.

#Libraries
library(tools)
library(ggplot2)
library(optparse)

#Command-line options
option_list <- list(
make_option(c("-s", "--sample"), help="Name of sequenced library"),
make_option(c("-r", "--reference"), help="Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)"),
make_option(c("-d", "--directory"), help="Directory (e.g., /projects/home/agombolay3/data/repository/Ribose-Map)")
)

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

for(i in opt$sample) {
        
        #Specify output directory and file
        directory <- file.path(opt$directory, "Results", opt$reference, opt$sample, "Coverage")
	input_files <- list.files(path=directory, pattern=".bed", full.names=TRUE, recursive=FALSE)
        
	for(file in input_files){
		
		#Plot only if files exist
        	if (file.exists(file)) {
            
        		#Specify dataset
            		data=read.table(file, sep="\t", header=FALSE)
			
			#Transform values on negative strand
			values <- ifelse(data$V6=='-',data$V3*-1,data$V3)

#############################################################################################################################
        		#Plot hotspots
        		myplot <- ggplot(data, aes(x=data[,2], y=values)) +

			geom_bar(stat = "identity",colour="black", fill="black")
			
			#Replace default theme
                	theme(panel.grid=element_blank(), panel.background=element_blank(),
			      axis.line=element_line(colour="black")) +
				
        		#Decrease space between plot and axes
        		scale_y_continuous(expand=c(0.015,0)) + scale_x_continuous(expand=c(0.015,0)) +
			
			#Add axes titles and specify font size
        		xlab("Chromosome Position") + ylab("rNMP Frequency") + theme(text=element_text(size=15))

#############################################################################################################################
#Save plot as PNG file
ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=20)

}
}
}
