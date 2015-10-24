#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-d", "--datasetFile"), help="input dataset description file"),
	make_option(c("-f", "--outFile"), help="output file having differentially bound regions"),
	make_option(c("-s", "--sessionFile"), help="output session file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$datasetFile) | is.null(opt$sessionFile) | is.null(opt$outFile)) {
	cat("\nProgram: diffPeaks.R (R script to identify differentially bound regions)\n")
	cat("Author: University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(session))

## load dataset description file, exit if not found
if(!file.exists(opt$datasetFile)) { cat("Cannot locate dataset description file\n"); q('no'); }
dataset <- dba(sampleSheet=opt$datasetFile, peakCaller="macs", scoreCol=8)

## calculate differentially bound regions (peaks)
data <- dba.count(dataset, minOverlap=2)
data <- dba.contrast(data, categories=DBA_CONDITION, minMembers=2)
data <- dba.analyze(data)
data.DB <- dba.report(data)
data.DB <- as.data.frame(data.DB)
write.table(data.DB, file=opt$outFile, sep="\t", quote=F, col.names=T, row.name=F)
write.table(data.DB[which(data.DB$Fold > 0),], file=sprintf("%s.onlywt", opt$outFile), sep="\t", quote=F, col.names=T, row.name=F)
write.table(data.DB[which(data.DB$Fold < 0),], file=sprintf("%s.onlyko", opt$outFile), sep="\t", quote=F, col.names=T, row.name=F)

## plot quality features
pdf("ma_plot.pdf")
dba.plotMA(data)
dev.off()

pdf("pca_plot.pdf")
par(mfrow=c(1,2))
dba.plotPCA(data, DBA_REPLICATE, attributes=DBA_CONDITION)
dba.plotPCA(data, contrast=1, th=0.05, attributes=DBA_CONDITION)
dev.off()

pdf("box_plot.pdf")
pvals <- dba.plotBox(data)
dev.off()

pdf("heatmap_plot.pdf")
corvals = dba.plotHeatmap(data, contrast=1, correlations=FALSE)
dev.off()

## determine and plot overlap rate between peaks from the input samples 
pdf("overlap_rate.pdf")
olap.rate <- dba.overlap(dataset, mode=DBA_OLAP_RATE)
plot(olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')
dev.off()

## determine peaks unique to the wt and ko
pdf("unique_count.pdf")
dba.overlap(dataset, dataset$masks$cebpa_wt, mode=DBA_OLAP_RATE)
dba.overlap(dataset, dataset$masks$cebpa_ko, mode=DBA_OLAP_RATE)
data <- dba.peakset(dataset, consensus=DBA_CONDITION, minOverlap=1)
dba.plotVenn(data, data$masks$Consensus)
dev.off()

data.OL <- dba.overlap(data, data$masks$Consensus)
write.table(as.data.frame(data.OL$onlyA), file=sprintf("%s.onlywtOccupancy", opt$outFile), sep="\t", quote=F, col.names=T, row.name=F)
write.table(as.data.frame(data.OL$onlyB), file=sprintf("%s.onlykoOccupancy", opt$outFile), sep="\t", quote=F, col.names=T, row.name=F)

## save the session for further use
save.session(opt$sessionFile)
