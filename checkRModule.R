#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"), help="package name"),
    make_option(c("-l", "--list"), help="instead list all installed packages in this output file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$input)) & (is.null(opt$list))) {
    cat("\nProgram: checkRModule.R (check if a R module is installed)\n")
    cat("Author: BRIC, University of Copenhagen, Denmark\n")
    cat("Version: 1.0\n")
    cat("Contact: pundhir@binf.ku.dk\n");
    print_help(parser)
    q()
}

if((!is.null(opt$input))) {
    is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

    cat(is.installed(opt$input))
} else {
    write.table(installed.packages()[,c(1,2,3,16)], file=opt$list, col.names=T, row.names=F, quote=F)
}
