#!/usr/bin/env Rscript

#© 2016 Alli Gombolay
#Author: Alli Lauren Gombolay
#E-mail: alli.gombolay@gatech.edu
#This program calculates the Chi Square statistic and its associated p-value.
#If the observed or expected counts are < 5, then categories will be collapsed.

#Load ggplot2
library("ggplot2")

#Argument 1 = dataset file
dataset <- commandArgs(trailingOnly = TRUE)[1]

#Argument 2 = PDF filename
filename <- commandArgs(trailingOnly = TRUE)[2]

#OBSERVED DATA

#Assign rNMP counts to variable
counts <- read.csv(file = dataset, sep = "\t")[ ,c('Counts')]

#Assign # of windows to variable
observed.windows <- read.csv(file = dataset, sep = "\t")[ ,c('Windows')]

#Sum # of windows in observed dataset
observed.total <- sum(observed.windows)

#Calculate proportions (# of windows/total)
observed.frequencies <- observed.windows/observed.total

#EXPECTED DATA

#Calculate Poisson distribution lambda value
lambda <- sum(observed.windows*counts/observed.total)

#Calculate frequencies based on lambda value
expected.frequencies <- c(dpois(counts, lambda))

#Calculate expected # of windows based on frequencies
expected.windows <- c(dpois(counts, lambda))*observed.total

#Determine if data need to be collapsed into few categories

#Rows with observed windows >=5
idx <- which(expected.windows >= 5)

#Row at which to start collapsing data
row.start <- length(expected.windows[idx])

#Number of rows of data in dataset
row.end <- length(readLines(dataset))-1

#Sum # of observed windows that need to be collapsed
collapsed.observed.sum <- sum(observed.windows[c(row.start:row.end)])

#Select rows of observed data that do not need to be collapsed
observed.data <- read.csv(file=dataset, sep = "\t")[ 1:row.start-1,]

#Append collapsed data to last row of dataset to create final observed dataset
observed.data.final <- rbind(observed.data, c(row.start-1, collapsed.observed.sum))
observed.windows.final <- observed.data.final[ ,c('Windows')]

#Sum # of expected windows that need to be collapsed
collapsed.expected.sum <- sum(expected.windows[c(row.start:row.end)])

#Select rows of expected data that do not need to be collapsed and afterward
#Append collapsed data to last row of dataset to create final expected dataset
expected.windows.final <- c(expected.windows[c(0:(row.start-1))],collapsed.expected.sum)

#CHI SQUARE TEST

#Calculate Chi Square Statistic from observed and expected number of windows in dataset
chi.square <- sum(((observed.windows.final - expected.windows.final)^2)/expected.windows.final)

#Calculate test statistic
#df=#categories-#estimated parameters-1
test.statistic <- qchisq(.95, df=row.start-1-1)

#Determine associated p-value
if (chi.square > test.statistic) {
  print("p-value < 0.05")
}

#Create dataframe
data <- data.frame(
  type = factor(c(replicate(row.start, "Poisson"), replicate(row.start, "Observed"))),
  counts = factor(c(seq(0,row.start-1,1),seq(0,row.start-1,1)), levels=seq(0,row.start-1,1)),
  windows = c(expected.windows.final,observed.windows.final)
)

#Plot observed and expected data
myplot <- ggplot(data=data, aes(x=counts, y=windows, fill=type))+geom_bar(stat="identity", position=position_dodge(width=0.8),
  width=0.6)+scale_fill_manual(values=c("#000000","#999999"))+xlab("rNMP count per 2.5 kb window")+ylab("Number of Genomic Windows")+
  theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+
  guides(fill=guide_legend(title=""))+scale_x_discrete(labels=seq(0,row.start-1,1))+scale_y_continuous(expand=c(0.015,0))+theme(text=element_text(size=14))

ggsave(filename=filename, plot=myplot)
