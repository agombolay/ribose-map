#!/usr/bin/env Rscript

#Copyright 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program plots rNMP and flanking nucleotide frequencies (100 bp upstream and downstream)
#rNMP frequencies = position 0, upstream frequencies = -100 -> -1, and downstream = +1 -> +100

#############################################################################################################################
#Load libraries
library(ggplot2); library(optparse)

#Command line options
option_list <- list(
	make_option(c("-s", "--sample"), help="Sequenced library name"),
	make_option(c("-d", "--directory"), help="Ribose-Map repository")")
)

#Get command line options, if -h invoked print help
opt <- parse_args(OptionParser(option_list=option_list))

#############################################################################################################################
#Specify output directory and file
path <- file.path(opt$directory, "Results", opt$sample, "Frequencies")
input <- list.files(path=path, pattern=".txt", full.names=T, recursive=F)

for(file in input){

	#Plot if file exists
	if (file.exists(input)) {
		
		#Plot regular and zoomed datasets
		for(i in c("Normal", "Zoomed")) {

#############################################################################################################################
			#Specify datasets to be used for each round of loop
			if (i=="Normal") {data=read.table(file, sep="\t", header=TRUE)}
			if (i=="Zoomed") {data=read.table(file, sep="\t", header=TRUE)[86:116,]}
    
			#Define variables to store nucleotide positions and frequency values
			position <- data$X; A <- data$A; C <- data$C; G <- data$G; T <- data$U.T

#############################################################################################################################
		#Create plot
		myplot <- ggplot(data, aes(x=position)) +

		#Remove and replace default background theme of plot
		theme(panel.grid=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"))+
			
		#Add axes titles and plot title
		xlab("Chromosome Position") + ylab("rNMP Frequency") +
		
		#Specify font size and remove grey background of legend key
		theme(text=element_text(size=20)) + theme(legend.key = element_blank()) +
		
		#Specify color values for each nucleotide and remove legend title
                scale_colour_manual(values=c("#CC79A7", "#56B4E9", "#E69F00", "#009E73"), name="") +
	
		#Plot data as scatterplot with connecting lines
                geom_line(aes(y=A, colour="A")) + geom_point(aes(y=A, colour="A")) + geom_line(aes(y=C, colour="C"))+ 
		geom_point(aes(y=C, colour="C")) + geom_line(aes(y=G, colour="G")) + geom_point(aes(y=G, colour="G"))+
                geom_line(aes(y=T, colour="U/T")) + geom_point(aes(y=T, colour="U/T"))

#############################################################################################################################
		#Save plot as PNG file
		ggsave(filename=file.path(output, paste(file_path_sans_ext(basename(file)), ".", i, ".png", sep="")), plot=myplot)
			
}
}
}
