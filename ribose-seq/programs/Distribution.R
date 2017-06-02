#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates the Chi Square statistic and its associated p-value.
#If the observed or expected counts are < 5, then categories will be collapsed.

#Load ggplot2
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
    for(j in c("all", "mito", "nucleus")) {
        
    #Specify output directory and file
    output <- file.path(opt$directory, "Ribose-Map", "Results", opt$reference, opt$sample, "Distribution")
    file <- file.path(output, paste(opt$sample, "-", "Counts", ".", opt$reference, ".", j, ".txt", sep=""))
    
    #Plot only if files exist
    if (file.exists(file)) {
#############################################################################################################################
        #OBSERVED DATA

        #Assign rNMP counts to variable
        number.rNMPs <- read.csv(file=file, sep="\t")[ ,c('rNMPs')]

        #Assign # of windows to variable
        observed.windows <- read.csv(file=file, sep="\t")[ ,c('Windows')]

#############################################################################################################################
        #EXPECTED DATA

        #Calculate Poisson distribution lambda value
        lambda <- sum(observed.windows*number.rNMPs)/sum(observed.windows)

        #Calculate expected # of windows based on frequencies
        expected.windows <- c(dpois(number.rNMPs, lambda))*sum(observed.windows)

#############################################################################################################################
        #COLLAPSE DATA IF EXPECTED COUNTS >=5

        #Rows with observed windows >=5
        idx <- which(expected.windows >= 5)

        #Print error message if idx equals zero
        if ((identical(idx, integer(0)))=="TRUE") {
            message("Plot cannot be created for", " ", i, " ", "(", j, ")")
        } else {

        #Row at which to start collapsing data
        row.start <- length(expected.windows[idx])

        #Number of rows of data in dataset
        row.end <- length(readLines(file))-1

        #Sum # of observed windows that need to be collapsed
        collapsed.observed.sum <- sum(observed.windows[c(row.start:row.end)])

        #Append collapsed data to last row of dataset to create re-binned expected dataset
        observed.windows.final <- c(observed.windows[1:row.start-1],collapsed.observed.sum)

        #Sum # of expected windows that need to be collapsed
        collapsed.expected.sum <- sum(expected.windows[c(row.start:row.end)])

        #Append collapsed data to last row of dataset to create re-binned expected dataset
        expected.windows.final <- c(expected.windows[1:row.start-1],collapsed.expected.sum)

#############################################################################################################################
        #CHI SQUARE GOODNESS OF FIT TEST FOR POISSON DISTRIBUTION

        #Calculate Chi Square Statistic from observed/expected number of windows in re-binned data
        chi.square <- sum(((observed.windows.final - expected.windows.final)^2)/expected.windows.final)

        #Calculate critical value at 95% confidence
        test.statistic <- qchisq(.95, df=row.start-1-1)

        #Determine if p-value is < 0.05
        if (chi.square > test.statistic) {
            message("p-value for", " ", i, " ", "(", j, ")", " ", "< 0.05")
        } else {
            message("p-value for", " ", i, " ", "(", j, ")", " ", ">= 0.05")
        }

#############################################################################################################################
        #PLOT OBSERVED VS. EXPECTED DISTRIBUTION OF RNMPS ACROSS GENOME

        #Create dataframe
        data <- data.frame(
            type = factor(c(replicate(row.start, "Poisson"), replicate(row.start, "Observed"))),
            rNMPs = factor(c(seq(0,row.start-1,1),seq(0,row.start-1,1)), levels=seq(0,row.start-1,1)),
            windows = c(expected.windows.final,observed.windows.final)
        )

        #Plot distributions
        myplot <- ggplot(data=data, aes(x=rNMPs, y=windows, fill=type)) +

        #Plot data as barcharts
        geom_bar(stat="identity",position=position_dodge(width=0.8), width=0.6) +

        #Remove and replace default background plot theme
        theme_bw() + theme(panel.border=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))+
    
        #Add axes titles and plot title (specified by user)
        xlab("rNMPs per Window") + ylab("Number of Windows") + ggtitle(opt$title) +
    
        #Specify font size and center title (if any) of plot on page
        theme(text=element_text(size=14)) + theme(plot.title=element_text(hjust=0.5)) +

        #Specify colors of bars, remove legend title, and decrease space between bars and x-axis
        scale_fill_manual(values=c("#000000", "#999999"), name="") + scale_y_continuous(expand=c(0.015,0))

#############################################################################################################################
#Save plot as PNG file
ggsave(filename=file.path(output, paste(opt$sample, "-", "Distribution", ".", j, ".png", sep="")), plot=myplot)

    }
    }
    }
message("Plotting of", " ", i , " ", "is complete")
}
