#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu

#1. Plots rNMP coverage at each chromosome position
#2. Saves plots as png files to appropriate directory

#Libraries
library(tools); library(ggplot2); library(optparse)

#Command-line options
option_list <- list(
	make_option(c("-s", "--sample"), help="Sequenced library name"),
	make_option(c("-d", "--directory"), help="Ribose-Map repository")
)

#Get command line options, if -h invoked print help
opt <- parse_args(OptionParser(option_list=option_list))

#############################################################################################################################
#Specify output directory and file
path <- file.path(opt$directory, "Results", opt$sample, "Coverage")
input <- list.files(path=path, pattern=".bed", full.names=T, recursive=F)
        
for(file in input){
		
	#Plot if file exists
        if (file.exists(file)) {
            
        	#Specify dataset
		data=read.table(file, sep="\t", header=FALSE)
			
		#Transform values on negative strand
		values <- ifelse(data$V6=='-',data$V3*-1,data$V3)

#############################################################################################################################
		#Create plot
		myplot <- ggplot(data, aes(x=data[,2], y=values))+
		
		#Create barchart and add titles for axes
        	geom_bar(stat="identity", fill="black") + xlab("Chromosome Position") + ylab("rNMP Frequency")+ 
		
		#Remove and replace default background theme of plot
		theme(panel.grid=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"))+
		
		#Decrease spacing between plot and axes and increase font
		theme(text=element_text(size=15)) + scale_y_discrete(expand=c(0.015,0)) + scale_x_discrete(expand=c(0.015,0))

#############################################################################################################################
#Save plot as PNG file
ggsave(filename=file.path(path, paste(file_path_sans_ext(basename(file)), ".png", sep="")), plot=myplot, width=20)

}
}
}
