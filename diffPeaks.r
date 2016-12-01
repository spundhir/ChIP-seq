#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-d", "--datasetFile"), help="input dataset description file"),
	make_option(c("-o", "--outDir"), help="output directory having differentially bound regions")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$datasetFile) | is.null(opt$outDir)) {
	cat("\nProgram: diffPeaks.R (R script to identify differentially bound regions)\n")
	cat("Author: University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
    cat("Format (description description file):\n");
    cat("SampleId,Tissue,Factor,Condition,Replicate,bamReads,Peaks\n");
    cat("suv39h1_control_Rep1,gmp,suv39h1,control,1,BOPNEQSSC1.bam,BOPNEQSSC1_peaks.broadPeak\n\n")
	q()
}

## load libraries
suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(session))

## load dataset description file, exit if not found
if(!file.exists(opt$datasetFile)) { cat("Cannot locate dataset description file\n"); q('no'); }

## create output directory, if does not exist
dir.create(file.path(opt$outDir), showWarnings = FALSE)

dataset <- dba(sampleSheet=opt$datasetFile, peakCaller="macs", scoreCol=8)

## calculate differentially bound regions (peaks)
data <- dba.count(dataset, minOverlap=1)
data <- dba.contrast(data, categories=DBA_CONDITION, minMembers=2)
data <- dba.analyze(data)

## report for all regions compared for differential binding analysis
data.DB <- dba.report(data, bCalled=T, th=1)
data.DB <- as.data.frame(data.DB)
write.table(data.DB, file=sprintf("%s/peaks.de", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)
write.table(data.DB[which(data.DB$Fold > 0),], file=sprintf("%s/peaks.de.onlywt", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)
write.table(data.DB[which(data.DB$Fold < 0),], file=sprintf("%s/peaks.de.onlyko", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)

## report for only differentially binding regions
data.DB <- dba.report(data, bCalled=T, th=0.1)
data.DB <- as.data.frame(data.DB)
write.table(data.DB, file=sprintf("%s/peaks.de.sig", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)
write.table(data.DB[which(data.DB$Fold > 0),], file=sprintf("%s/peaks.de.onlywt.sig", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)
write.table(data.DB[which(data.DB$Fold < 0),], file=sprintf("%s/peaks.de.onlyko.sig", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)

## plot quality features
pdf(sprintf("%s/ma_plot.pdf", opt$outDir))
dba.plotMA(data)
dev.off()

pdf(sprintf("%s/pca_plot.pdf", opt$outDir))
par(mfrow=c(1,2))
dba.plotPCA(data, DBA_REPLICATE, attributes=DBA_CONDITION)
dba.plotPCA(data, contrast=1, th=0.05, attributes=DBA_CONDITION)
dev.off()

pdf(sprintf("%s/box_plot.pdf", opt$outDir))
pvals <- dba.plotBox(data)
dev.off()

pdf(sprintf("%s/heatmap_plot.pdf", opt$outDir))
corvals = dba.plotHeatmap(data, contrast=1, correlations=FALSE)
dev.off()

## determine and plot overlap rate between peaks from the input samples 
pdf(sprintf("%s/overlap_rate.pdf", opt$outDir))
olap.rate <- dba.overlap(dataset, mode=DBA_OLAP_RATE)
plot(olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets')
dev.off()

## determine peaks unique to the wt and ko
## NOTE: above differential binding analysis is solely based on tag count between two conditions
##       overlap analysis below is based on if peaks are called by peak caller (macs2) between two conditions
pdf(sprintf("%s/unique_count.pdf", opt$outDir))
dba.overlap(dataset, dataset$masks$cebpa_wt, mode=DBA_OLAP_RATE)
dba.overlap(dataset, dataset$masks$cebpa_ko, mode=DBA_OLAP_RATE)
data <- dba.peakset(dataset, consensus=DBA_CONDITION, minOverlap=1)
dba.plotVenn(data, data$masks$Consensus)
dev.off()

data.OL <- dba.overlap(data, data$masks$Consensus)
write.table(as.data.frame(data.OL$onlyA), file=sprintf("%s/peaks.de.onlywtOccupancy", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)
write.table(as.data.frame(data.OL$onlyB), file=sprintf("%s/peaks.de.onlykoOccupancy", opt$outDir), sep="\t", quote=F, col.names=T, row.name=F)

## save the session for further use
save.session(sprintf("%s/peaks.de.session", opt$outDir))
