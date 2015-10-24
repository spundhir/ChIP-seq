#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing results from bed2superEnhancer script"),
	make_option(c("-o", "--outFile"), help="output file containing predictions"),
	make_option(c("-t", "--threshold"), help="total expression threshold at x-axis above which to consider an enhancer as super enhancer (computed by bed2superEnhancer.py)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outFile)) {
	cat("\nProgram: bed2direction.R (R script to predict super enhancers)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(ggplot2))

## function to scale values between 0 and 1
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

data <- read.table(opt$inFile)
data <- data[order(data[,9]),]
data$length <- data$V3-data$V2
data$superEnhancer <- "No"
#data$normExprScaled <- scale01(data$V4)
#data$rawExprScaled <- scale01(data$V5)
data[as.integer(opt$threshold):length(data[,1]),]$superEnhancer <- "Yes"
color <- c(rep("#000000", length(which(data$superEnhancer=="No"))), rep("red", length(which(data$superEnhancer=="Yes"))))
pdfFile <- paste(gsub("\\..*", "", opt$outFile), ".pdf", sep="")
pdf(pdfFile)
barplot(data$V8, ylab="Normalized ChIP-seq signal", xlab="Enhancers ranked by total ChIP-seq signal", cex=1.5, cex.axis=1.5, cex.lab=1.5, col=color, border=NA)
dev.off()

data <- data[order(-data[,9]),]
write.table(data[,c(1:11)], opt$outFile, sep="\t", row.names=F, col.names=F, quote=F)
