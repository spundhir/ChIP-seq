#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inputFile"), help="input file from homer"),
	make_option(c("-o", "--outputFile"), help="output pdf file"),
	make_option(c("-m", "--mode"), help="different modes available to plot the data")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inputFile) | is.null(opt$outputFile) | is.null(opt$mode)) {
	cat("\nProgram: plotHomerResults.r (R script to plot results from homer)\n")
	cat("Author: University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

if(opt$mode==1) {
    data <- read.table(opt$inputFile, sep="\t", header=T)
    colnames(data) <- gsub("_.*|\\..*|.Uniprobe.*|.*BestGuess.", "", colnames(data))
    pdf(opt$outputFile)
    par(mfrow=c(2,1))
    for(i in seq(2, length(colnames(data))-6, 3)) { barplot(data[,i], names.arg=data[,1], las=2, xlab="distance to peak", ylab="density", main=colnames(data)[i]); }
    dev.off()
}
