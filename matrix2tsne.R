#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing data matrix (Col 1: class followed by values"),
    make_option(c("-o", "--outFile"), help="output pdf file"),
    make_option(c("-s", "--sessionFile"), help="output session file"),
    make_option(c("-p", "--perplexity"), default=30, help="perplexity value numeric (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outFile)) {
	cat("\nProgram: matrix2tsne.R (R script to plot TSNE from matrix)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(Rtsne))

ggpca <- function(data, groups = NULL, main = "", palette = "Paired", ...){

    # Groups flag
    gFlag <- TRUE

    if (is.null(groups)){
        groups <- rep("sample", nrow(pca$x))
        gFlag <- FALSE
    }

    data <- data.frame(groups, data)

    colnames(data) <- c("groups", "x", "y")

    #Prevent R CMD check "no global definition.." NOTE
    x <- y <- NULL

    if (gFlag) {
        p <- ggplot(data, aes(x = x, y = y, label = groups,
            colour = factor(groups)))
            p <- p + scale_color_brewer(palette = palette) + labs(colour='Groups')
    } else {
        p <- ggplot(data, aes(x = x, y = y))
    }

    p + geom_point()
}

if(opt$inFile=="stdin") {
    data <- read.table(stdin(), sep="\t")
} else {
    data <- read.table(opt$inFile, sep="\t")
}
#data <- unique(data)
#set.seed(1)
#tsne_out <- Rtsne(as.matrix(log(data[,c(2:ncol(data))])), check_duplicates=F, perplexity=as.numeric(opt$perplexity))
#ggpca(tsne_out$Y, groups=data[,1]) + theme_bw()
test <- unique(data[,c(5,16:19)])
tsne_out <- Rtsne(as.matrix(log(test[,c(2:5)]))) # Run TSNE
ggpca(tsne_out$Y, groups = data[,5]) + theme_bw()
ggsave(dpi=300, width=10, height=10, filename=opt$outFile, useDingbats=FALSE)
if(!is.null(opt$sessionFile)) {
    save.session(opt$sessionFile)
}
q()
