#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--queryDF"), help="query data frame as file whose distribution of row values are to be compared"),
	make_option(c("-j", "--subjectDF"), help="subject data frame as file against which distribution of row values are to be compared ")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$queryDF) | is.null(opt$subjectDF)) {
	cat("\nProgram: ksTestOnDF.R (R script to perform Kolmogorov-Smirnov test iteratively between all rows of two data frames)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))

f <- function(x, y, ..., alternative = c("two.sided", "less", "greater"), exact = NULL){ 
    #w <- getOption("warn") 
    #options(warn = -1)  # ignore warnings 
    p <- ks.test(x, y, ..., alternative = alternative, exact = exact)$p.value 
    #options(warn = w) 
    p 
}

scale01 <- function(x){
    (x-min(x))/(max(x)-min(x))
}

data1 <- read.table(opt$queryDF)
data2 <- read.table(opt$subjectDF)
data1$V11 <- data1$V7-data1$V8
data1$V12 <- scale01(data1$V10)-scale01(data1$V9)
data2$V11 <- data2$V7-data2$V8
data2$V12 <- scale01(data2$V10)-scale01(data2$V9)

ss <- apply(data1[,c(7:12)], 1, function(y) mean(apply(data2[which(data2$V4=="SS"),c(7:12)], 1, function(x) f(x,y))))
su <- apply(data1[,c(7:12)], 1, function(y) mean(apply(data2[which(data2$V4=="SU"),c(7:12)], 1, function(x) f(x,y))))
us <- apply(data1[,c(7:12)], 1, function(y) mean(apply(data2[which(data2$V4=="US"),c(7:12)], 1, function(x) f(x,y))))
uu <- apply(data1[,c(7:12)], 1, function(y) mean(apply(data2[which(data2$V4=="UU"),c(7:12)], 1, function(x) f(x,y))))

all <- data.frame(ss, su, us, uu)
write.table(all, "test.ks", sep="\t", quote=F, row.names=F, col.names=F)
save.session("test.ks.session")
