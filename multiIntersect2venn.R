#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input BED file created using multiIntersectBed (can be stdin)"),
	make_option(c("-o", "--outFile"), help="output pdf file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outFile)) {
	cat("\nProgram: multiIntersect2venn.R (R script to plot venn diagram from multiIntersectBed results)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library("RColorBrewer"))

if(opt$inFile=="stdin") {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$inFile)
}
vec <- as.vector(unique(data[!grepl(",", data$V5),]$V5))
data$id <- sprintf("%d_%d", data$V2, data$V3)
lst <- list()
k <- 1
for(i in vec) { 
    l <- data[grep(i, data$V5),]$id
    lst[[k]] <- l
    k=k+1
}
names(lst) <- vec
col <- brewer.pal(length(vec)+1, "Spectral")
col <- col[1:length(vec)]
venn.plot <- venn.diagram(lst, fill=col, NULL)
pdf(opt$outFile)
grid.draw(venn.plot)
dev.off()
