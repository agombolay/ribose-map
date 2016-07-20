# nuc.freqs.R
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 559 $'
#
# modmap pipeline for plotting nucleotide frequencies 

library(ggplot2)
library(RColorBrewer)
library(Cairo)
library(optparse)

option_list <- list(
    make_option(c("-n", "--name"), action="store", 
                default="", help="sample name"),
    make_option(c("-d", "--directory"), action="store",
                default=".", help="output directory [default %default]"),
    make_option(c("--offsetmax"), action="store",
                default="15", help="offset max [default %default]"))

parser <- OptionParser(usage="%prog [options] infile",
                       option_list=option_list)

arguments <- parse_args(parser, positional_arguments = 1)

opt <- arguments$options
infile <- arguments$args

sample.name <- opt$name
output.dir <- opt$directory

offset.max <- as.integer(opt$offsetmax)
offset.min <- as.integer(offset.max) * -1

COLNAMES <- c('nuc','offset','region.size','count',
              'freq','type','total.sites')
df <- read.table(infile, col.names=COLNAMES)

if (nrow(df) == 0) {
    warning("empty data frame")
    quit(status=0)
}
head(df)

# color interpolation to expand palette
# http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot.nuc.freq <- function(df, cur.size, offset.min, offset.max, ... ) {

    # subset the data
    df <- subset(df, region.size == cur.size)

    # total.sites should all be the same, just take the first
    num.sites <- df$total.sites[1]

    sum.counts <- sum(df$count)
  
    # number of colors needed for the palette
    num.colors <- length(unique(df$nuc))

    gp <- ggplot(data = df,
                     aes(nuc, freq, offset))

    gp <- gp + geom_line(aes(color=factor(nuc),
                             x = offset, y = freq))

    gp <- gp + geom_point(aes(x = offset, y = freq, 
                              color = nuc, size = 3),
                          show_guide = FALSE)

    gp <- gp + facet_wrap(~type, scales = "free_y")
  
    gp <- gp + theme_bw()
    gp <- gp + theme(legend.position = 'bottom')

    if (num.colors > 4) {
        gp <- gp + scale_color_manual(values = getPalette(num.colors))
        gp <- gp + guides(color = guide_legend(title = NULL,
                                               nrow = 4))
    } else {
        gp <- gp + scale_color_brewer(palette="Set1")
        gp <- gp + guides(color = guide_legend(title = NULL,
                                               nrow = 1))
    }

    gp <- gp + xlim(offset.min, offset.max)

    # axis labels 
    gp <- gp + xlab('Position')
    gp <- gp + ylab('Frequency')

    # add title
    title.top = paste('modmap nucleotide-frequency\n(sample ',
                      sample.name, ' region size ', cur.size, ')', sep='')
    title.bottom = paste('n.sites = ', num.sites, ' n.counts = ',
                         sum.counts, sep='')

    title = paste(title.top, title.bottom, sep='\n')
    gp <- gp + ggtitle(title)

    return(gp)
}

uniq.sizes = unique(df$region.size)

for (idx in 1:length(uniq.sizes)) {

    cur.size <- uniq.sizes[idx]

    gp.nuc.freq <- ggplot.nuc.freq(df, cur.size, offset.min, offset.max)

    # write the file
    pdf.filename <- paste(output.dir, '/', 'modmap.nuc.freq.region.',
                          cur.size, '.', 
                          "min.", offset.min, ".", 
                          "max.", offset.max,
                          ".", sample.name, '.pdf', sep='')

    ggsave(filename = pdf.filename, plot = gp.nuc.freq,
           height = 8.5, width = 11,
           device = CairoPDF)
}
