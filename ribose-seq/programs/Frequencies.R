#!/usr/bin/env Rscript

#Copyright 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots rNMP and flanking nucleotide frequencies (100 bp upstream and downstream)
#rNMP frequencies = position 0, upstream frequencies = -100 -> -1, and downstream = +1 -> +100

#############################################################################################################################
#Load libraries
library(optparse)
library(ggplot2)

#Command line options
option_list <- list(
make_option(c("-s", "--sample"), help="Sample name(s) (e.g., FS1, FS2, FS3)"),
make_option(c("-r", "--reference"), help="Reference genome (e.g., sacCer2, pombe, ecoli, mm9, hg38)"),
make_option(c("-d", "--directory"), help="Directory (e.g., /projects/home/agombolay3/data/repository/Ribose-Map)")
)

#Get command line options, if -h encountered print help
opt <- parse_args(OptionParser(option_list=option_list))

#############################################################################################################################
for(i in opt$sample) {
	for(j in c("mito", "nucleus")) {

		#Specify output directory and file
		output <- file.path(opt$directory, "Results", opt$reference, opt$sample, "Frequencies")
		file <- file.path(output, paste(opt$sample, "-", "Frequencies", ".", j, ".txt", sep=""))

 		#Plot only if files exist
        	if (file.exists(file)) {
			#Plot regular and zoomed datasets
                	for(k in c("Regular", "Zoomed")) {

#############################################################################################################################
                    	#Specify datasets to be used for each round of loop
                    	if (k=="Regular") {data=read.table(file, sep="\t", header=TRUE)}
                    	if (k=="Zoomed") {data=read.table(file, sep="\t", header=TRUE)[86:116,]}
    
                    	#Define variables to store nucleotide positions and frequency values
                    	position <- data$X; A <- data$A; C <- data$C; G <- data$G; T <- data$U.T

#############################################################################################################################
                    	#Plot frequencies
                    	myplot <- ggplot(data=data, aes(x=position)) +
    
                    	#Replace default theme
        		theme(panel.grid=element_blank(),
              		panel.background=element_blank(),
			axis.line=element_line(colour="black")) +
			
			#Add axes titles and plot title
			xlab("Position") + ylab("Frequency") +
				
			#Plot data as scatterplot with connecting lines
                    	geom_line(aes(y=A, colour="A")) + geom_point(aes(y=A, colour="A")) +
                    	geom_line(aes(y=C, colour="C")) + geom_point(aes(y=C, colour="C")) +
                    	geom_line(aes(y=G, colour="G")) + geom_point(aes(y=G, colour="G")) +
                    	geom_line(aes(y=T, colour="U/T")) + geom_point(aes(y=T, colour="U/T")) +
				
                    	#Specify font size and center title (if any) of plot on page
                    	theme(text=element_text(size=20)) + theme(legend.key = element_blank()) +

                    	#Specify color values for each nucleotide and remove legend title
                    	scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="")
														     
#############################################################################################################################
#Save plot as PNG file
ggsave(filename=file.path(output, paste(opt$sample, "-", "Frequencies", "-", k, ".", j, ".png", sep="")), plot=myplot)
			
}
}
}
}
