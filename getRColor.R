#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--count"), help="number of colors required"),
  make_option(c("-c", "--class"), default="RdYlBu", help="color brewer class (default: %default)"),
  make_option(c("-p", "--paired"), help="two consecutive colors of same shade", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$count)) {
  cat("\nProgram: getRColor.R (get specified number of colors in hexadecimal)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

## define color class, if not defined or if paired is set
if(!is.null(opt$paired)) {
    opt$class <- "Paired"
}

## create color codes
if(as.numeric(opt$count) <= 11) {
    cat(brewer.pal(as.numeric(opt$count), opt$class))
} else {
    count <- 10
    pal <- colorRampPalette(brewer.pal(count, opt$class))
    cat(pal(opt$count))
}
#cat("\n")
q()        
