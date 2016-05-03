#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file created by nfrDynaAna script"),
    make_option(c("-n", "--minFreq"), default="50", help="minimum frequency of NFRs in each nfr dynamic class (defaut=%default)"),
    make_option(c("-d", "--diffFreq"), default="0", help="minimum difference in enrichment between categories (defaut=%default)"),
	make_option(c("-o", "--outPdfFile"), help="output pdf image file"),
    make_option(c("-f", "--onlyDiffFreq"), action="store_true", help="less stringent criteria. use only difference in frequency"),
    make_option(c("-v", "--plotSim"), action="store_true", help="instead plot motifs similar in frequency"),
    make_option(c("-r", "--rescale"), action="store_true", help="rescale color bar between min and max value")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outPdfFile)) {
	cat("\nProgram: motifDynAna.R (R script to plot motif dynamics)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(session))

data <- read.table(opt$inFile)
data$V13 <- log2((data$V8+0.01)/(data$V10+0.01))
data$V14 <- gsub("^.*BestGuess:", "", data$V2)
data$V14 <- gsub("\\(.*", "", data$V14)
data$V15 <- gsub("^.*BestGuess:", "", data$V2)
data$V15 <- as.numeric(gsub("\\)", "", gsub("^.*\\(", "", data$V15)))
data <- data[which(data$V11>=as.numeric(opt$minFreq)),]
no_rows=nrow(data)/length(unique(data$V1))
tf_info <- as.data.frame(data[1:no_rows, c(14,15)])
dat <- matrix(data$V8, nrow=no_rows)
pat <- matrix(data$V4, nrow=no_rows)
## use log fold change
mat <- matrix(data$V13, nrow=no_rows)
## use q-value
#mat <- matrix(data$V6, nrow=no_rows)
data$V1 <- sprintf("%s_%s", data$V1, data$V11)
colnames(mat) <- as.vector(unique(data$V1))
row.names(mat) <- data[1:no_rows,2]
if(!is.null(opt$onlyDiffFreq)) {
    #myCol <- brewer.pal(9, "OrRd")
    myCol <- rev(brewer.pal(11, "RdBu"))
    sig_rows <- which(apply(mat, 1, function(x) max(x)-min(x)>as.numeric(opt$diffFreq)))
} else if(!is.null(opt$plotSim)) {
    sig_rows <- which(apply(mat, 1, function(x) max(x) > 0.1 & min(x) > 0.1 & max(x)-min(x) <= as.numeric(opt$diffFreq)))
    myCol <- (brewer.pal(9, "OrRd"))
} else {
    myCol <- rev(brewer.pal(11, "RdBu"))
    sig_rows <- which(apply(mat, 1, function(x) max(x) > 0 & min(x) < 0 & max(x)-min(x) > as.numeric(opt$diffFreq)))
    #sig_rows <- which(apply(mat, 1, function(x) max(x) > as.numeric(opt$diffFreq) & min(x) < -1*as.numeric(opt$diffFreq) & max(x)-min(x) > 1))
}
#sig_rows <- which(apply(dat, 1, function(x) max(x)>3))

if(length(sig_rows)>2) {
    pdf(opt$outPdfFile, height=15)
    #breaks <- as.vector(summary(as.vector(mat[sig_rows,])))
    #len=2
    #breaks1 <- seq(breaks[1], breaks[2], length=len)
    #breaks2 <- seq(breaks[2], breaks[3], length=len)
    #breaks3 <- seq(breaks[3], breaks[4], length=len)
    #breaks4 <- seq(breaks[4], breaks[5], length=len)
    #breaks5 <- seq(breaks[5], breaks[6], length=len)
    #breaks <- c(breaks1[1:length(breaks1)], breaks2[2:length(breaks2)],
    #            breaks3[2:length(breaks3)], breaks4[2:length(breaks4)],
    #            breaks5[2:length(breaks5)])
    if(!is.null(opt$rescale)) {
        breaks <- seq(min(mat[sig_rows,]), max(mat[sig_rows,]), by=0.25)
        #myCol <- colorpanel(n=length(breaks)-1,low="blue",mid="white",high="red")
        myCol <- rev(brewer.pal(length(breaks)-1, "RdBu"))
        heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(15,25), cexCol=1, cexRow=1, breaks=breaks)
    } else {
        heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(15,25), cexCol=1, cexRow=1)
        #heatmap.2(mat, trace="none", col=myCol, margins=c(15,20), cexCol=1, cexRow=1)
    }
    dev.off()
} else if(length(sig_rows)>=1) {
    cat("\nNumber of significant motifs are too low for a heatmap. Instead writing results to motif_dynamics.txt\n\n")
    if(length(sig_rows)==1) {
        write(rownames(mat)[sig_rows], file="motif_dynamics.txt", sep="\t")
        write.table(mat[sig_rows,], "motif_dynamics.txt", quote=F, append=TRUE, sep="\t")
    } else {
        write.table(mat[sig_rows,], "motif_dynamics.txt", quote=F, append=FALSE, sep="\t")
    }
} else {
    cat("\nNo significant motif is found\n")
}
save.session("test.session")

## old code (v2.0)
#data <- read.table(opt$inFile)
#data$V10 <- apply(data, 1, function(x) binom.test(as.numeric(x[6]), as.numeric(x[7]), as.numeric(x[8])/as.numeric(x[9]), alternative="g")$p.value)
#no_rows=nrow(data)/length(unique(data$V5))
#dat <- matrix(data$V2, nrow=no_rows)
## use log fold change
#mat <- matrix(data$V4, nrow=no_rows)
## use difference between real and background percentages
#mat <- matrix(data$V2-data$V3, nrow=no_rows)
## use p-value computed using binomial test
#mat <- matrix(data$V10, nrow=no_rows)
#colnames(mat) <- as.vector(unique(data$V5))
#row.names(mat) <- data[1:no_rows,1]
#myCol <- brewer.pal(9, "Blues")
#sig_rows <- which(apply(dat, 1, function(x) max(x))>0.33)
#pdf(opt$outPdfFile)
#heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(12,20), cexCol=1, cexRow=1)
#dev.off()

## old code (v1.0)
#data$count <- 1
#dat <- dcast(data, V2 ~ V1)
#dat[is.na(dat)] <- 0
#cols <- ncol(dat)
#dat$mean <- apply(dat[,c(2:cols)],1, mean)
#dat <- dat[with(dat, order(-mean)),]
#mat <- data.matrix(dat)[,c(2:cols)]
#rownames(mat) <- dat$V2
#mycol <- colorpanel(n=99,low="white",high="black")
#colnames(mat) <- gsub("granulocytes", "grn", colnames(mat))
#pdf(opt$outPdfFile)
#heatmap.2(mat, trace="none", col=mycol, margins=c(11,10), cexCol=1)
#dev.off()
