#!/usr/bin/env Rscript

## help(distribution) to see other distributions
## input is a file containing even number of columns
## assumes consecutive two columns as values for each binomial test
## example: 20 30 40 12;
## binom.test(c(20,30), 0.5, alternative="g")$p.value
## binom.test(c(40,12), 0.5, alternative="g")$p.value

suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
    make_option(c("-i", "--countFile"), help="input read count file (format: <int> <int> (.. <int> <int>..) | can be stdin)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$countFile)) {
    cat("\nProgram: binom.R (R script to perform binomial test during pol2 pausing analysis)\n")
    cat("Author: BRIC, University of Copenhagen, Denmark\n")
    cat("Version: 1.0\n")
    cat("Contact: pundhir@binf.ku.dk\n");
    print_help(parser)
    q() 
}

if(identical(opt$countFile, "stdin")==T) {
    df <- read.table(file("stdin"), header=F)
} else {
    df <- read.table(opt$countFile, header=F)
}
j <- 1
for(i in seq(1,ncol(df),by=2)) {
    df <- cbind(df, p.adjust(apply(df, 1, function(x) binom.test(c(round(as.numeric(x[i])*100,0), round(as.numeric(x[i+1])*100,0)), 0.5, alternative="g")$p.value), "BH"))
    colnames(df)[ncol(df)] <- j
    j <- j +1
}
write.table(df, "", quote=F, col.names=F, row.names=F, sep="\t")
q()
