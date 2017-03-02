#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inFile"), help="input file containing values to fit in matrix format (can be stdin)"),
  make_option(c("-d", "--distribution"), default="pois", help="name of distribution to fit (pois, nbinom, gamma, weibull, nnorm, gumbel, binom, geom, hyper) (if multiple separate them by a comma) (default=%default)"),
  make_option(c("-c", "--discrete"), default=T, help="if distribution of discrete or continuous (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: matrix2fitDistr.R (R script to fit distributions)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(fitdistrplus))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ismev))
suppressPackageStartupMessages(library(session))

if(identical(opt$inFile, "stdin")==T) {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$inFile)
}

fitDistr <- function(values,distr) {
  ## extract input distribution names
  distr=as.vector(unlist(strsplit(distr, ",")))

  if(length(values)>=10) { 
    ## fit each input distribution to the input numeric values
    model <- lapply(distr, function(x) { fitdist(as.integer(values), x, method="mme") })

    for (i in 1:length(model)) {
        if(distr[i]=="pois") {
            param <- as.vector(coef(model[[i]]))
            pvalue <- laply(as.integer(values), function(x) { ppois(x, lambda = param[1], lower.tail = F) })
        }
        else if(distr[i]=="nbinom") {
            param <- as.vector(coef(model[[i]]))
            pvalue <- laply(as.integer(values), function(x) { pnbinom(x, size = param[1], mu=param[2], lower.tail = F) })
        }
        else if(distr[i]=="geom") {
            param <- as.vector(coef(model[[i]]))
            pvalue <- laply(as.integer(values), function(x) { pgeom(x, param[1], lower.tail = F) })
        }
    }
  } else {
    pvalue <- 1
  }
  return(sprintf("%E", pvalue))
}

# data <- read.table("nfr_dynamics_count/matrix/NFR_DYNAMICS_SIG.stat")
# distr <- "nbinom"
# descdist(as.numeric(data[,22]))
# model <- fitdist(as.numeric(data[,c(25)]), distr, method="mme")
# plot(model)
# cdfcomp(model, legendtext=distr)
# qqcomp(model, legendtext=distr)
# model$aic
# data.pVal <- cbind(data, apply(data[,c(10:25)], 2, function(x) fitDistr(x, "pois")))

## fit the distribution
data.pVal <- apply(data, 2, function(x) fitDistr(x, opt$distribution))

## write output containing pvalue information
write.table(data.pVal, "", sep="\t", row.names = F, col.names = F, quote = F)

q()