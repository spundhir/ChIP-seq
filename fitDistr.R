#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing numeric value for fit (first column) (can be stdin)"),
    make_option(c("-o", "--outDir"), help="output directory"),
    make_option(c("-d", "--distribution"), default="pois", help="name of distribution to fit (pois, nbinom, gamma, weibull, nnorm, gumbel, binom, geom, hyper) (if multiple separate them by a comma) (default=%default)"),
    make_option(c("-c", "--discrete"), default=T, help="if distribution of discrete or continuous (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outDir)) {
	cat("\nProgram: fitDistr.R (R script to fit distributions)\n")
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

## create output directory, if does not exist
dir.create(file.path(opt$outDir), showWarnings = FALSE)

## extract input distribution names
distr=as.vector(unlist(strsplit(opt$distr, ",")))

if(nrow(data)>=10) { 
    ## analyze each model fit
    outPdf <- sprintf("%s/plots.pdf", opt$outDir)
    pdf(outPdf)
    
    ## check as to which distribution is the best
    descdist(as.integer(data[,1]))

    ## fit each input distribution to the input numeric values
    model <- lapply(distr, function(x) { fitdist(as.integer(data[,1]), x, method="mme") })

    for (i in 1:length(model)) {
        if(distr[i]=="pois") {
            param <- as.vector(coef(model[[i]]))
            plot(model[[i]])
            data$pois <- laply(as.integer(data[,1]), function(x) { ppois(x, lambda = param[1], lower.tail = F) })

            ## write AIC value
            outFile <- sprintf("%s/AIC", opt$outDir)
            write(sprintf("POIS: %s", model[[i]]$aic), outFile, append=T)
        }
        else if(distr[i]=="nbinom") {
            param <- as.vector(coef(model[[i]]))
            plot(model[[i]])
            data$pnorm <- laply(as.integer(data[,1]), function(x) { pnbinom(x, size = param[1], mu=param[2], lower.tail = F) })

            ## write AIC value
            outFile <- sprintf("%s/AIC", opt$outDir)
            write(sprintf("NBINOM: %s", model[[i]]$aic), outFile, append=T)
        }
        else if(distr[i]=="geom") {
            param <- as.vector(coef(model[[i]]))
            plot(model[[i]])
            data$geom <- laply(as.integer(data[,1]), function(x) { pgeom(x, param[1], lower.tail = F) })

            ## write AIC value
            outFile <- sprintf("%s/AIC", opt$outDir)
            write(sprintf("GEOM: %s", model[[i]]$aic), outFile, append=T)
        }
    }

    #gofstat(model[[i]], fitnames=distr)
    cdfcomp(model, legendtext=distr)
    qqcomp(model, legendtext=distr)
    dev.off()

    ## HOW TO EVALUATE WHICH MODEL IS BETTER
    # one having lower AIC value (model[[i]]$aic
    # more details: http://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
} else {
    data$pois <- 1
}
 
## write output file containing pvalue information
outTable <- sprintf("%s/table", opt$outDir)
write.table(data[order(-data[,1]),], outTable, sep="\t", row.names = F, col.names = F, quote = F)

## write output session file
outSession <- sprintf("%s/session", opt$outDir)
save.session(outSession)
q()

## MISCELLANEOUS INFO ON DISTRIBUTION FITTING
library(ismev)
t <- gum.fit(y)
gum.diag(t)

t <-fitdistr(as.integer(y), densfun="negative binomial")
pnbinom(50, size = 0.18, mu=14.99, lower.tail = F)

t <-fitdistr(as.integer(y), densfun="poisson")
ppois(20, lambda=as.numeric(as.vector(t)[1]$estimate[1]), lower.tail = F)

set.seed(1)
x = seq(-2, 8, .01)
y = rnbinom(length(x), mu=exp(x), size=1)
fit = glm.nb(y ~ 1)
pnbinom(0.00001, size=fit$theta, prob=0.05)

predicted.y = predict(fit, newdata=data.frame(x=5), type="response")
dnbinom(100, mu=5000, size=fit$theta)
prob = function(newx, newy, fit) {
    dnbinom(newy, mu=predict(fit, newdata=data.frame(x=newx), type="response"), size=fit$theta)
}
