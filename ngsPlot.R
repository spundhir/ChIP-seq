#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--ngsPlotSession"), help="input session file corresponding to previous ngsPlot run"),
  make_option(c("-l", "--logScale"), help="plot histone profile in log scale", action="store_true"),
  make_option(c("-o", "--outfile"), help="output pdf file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$ngsPlotSession) | is.null(opt$outfile))) {
  cat("\nProgram: ngsPlot.R (plot scaled histone profiles computed using ngsPlot)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(caTools))
progpath <- "/home/pundhir/software/ngsplot-master"
source(file.path(progpath, 'lib', 'plotlib.r'))

load(opt$ngsPlotSession)
regcovMat <- as.matrix(runmean(regcovMat, k=mw, alg='C', endrule='mean'))
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
pdf(opt$outfile)
ymax <- 0
mw <- 0

if(is.null(opt$logScale)) {
    regcovMat <- apply(regcovMat, 2, function(x) scale01(x))
    plotmat(regcovMat, ctg.tbl$title, ctg.tbl$color, bam.pair, xticks, pts, m.pts, f.pts, pint, shade.alp, confiMat, mw, prof.misc)
} else {
    plotmat(log(regcovMat), ctg.tbl$title, ctg.tbl$color, bam.pair, xticks, pts, m.pts, f.pts, pint, shade.alp, confiMat, mw, prof.misc)
}
dev.off()

#regcovMat <- as.data.frame(regcovMat)
#regcovMat$id <- rownames(regcovMat)
#regcovMat$id <- factor(regcovMat$id, levels=c(regcovMat$id), ordered=TRUE)
#regcovMat.melt <- melt(regcovMat, id="id")
#ggplot(regcovMat.melt, aes(x=id, y=scale01(value), col=variable, group=variable)) + geom_line() + scale_color_manual(values=color)



## OLD CODE
#setwd("~/Lobster/project/chip-seq-analysis/analysis/peak_calling_histone/human/h3k4me1_helas3/ngsPlot")
#dirs = c("~/Lobster/project/chip-seq-analysis/analysis/peak_calling_histone/human/h3k4me1_helas3/ngsPlot")
#setwd(opt$ngsplot)
#dirs = c(opt$ngsplot)
#avgProfileAll <- data.frame()
#for(dir in dirs) {
#  avgProfile <- data.frame()
#  i <- 1
#  colnames <- list()
#  for(subdir in list.dirs(dir, recursive=F)) { 
#    load(paste(subdir, "combined/avgprof.RData", sep="/"))
#    if(i==1) {
#      avgProfile <- rbind(avgProfile, regcovMat)
#    } else {
#      avgProfile[,i] <- regcovMat[,1]
#    }
#    colnames = c(colnames, gsub("_.+_", "_", colnames(regcovMat)))
#    i <- i + 1
#  }
#  colnames(avgProfile) <- colnames
#  avgProfileAll <- rbind(avgProfileAll, avgProfile)
#}

#pdf(opt$outfile)
#library(gridExtra)
#library(ggplot2)
#colnames <- unlist(colnames)
#color <- c("#5F192C", "#C86368", "#EDAA4E", "#FCC059", "#4A789C", "#85C1F5")
#RPM: Read count Per Million mapped reads
#ggplot(data=avgProfileAll[1:101,], aes(c(1:101))) + xlab(NULL) + scale_x_discrete(breaks = seq(1, 101, by=25), labels=c("-1000","-500","0","500", "1000")) + geom_line(aes(y=get(colnames[1]), colour=colnames[1])) + geom_line(aes(y=get(colnames[2]), colour=colnames[2])) + geom_line(aes(y=get(colnames[3]), colour=colnames[3])) + geom_line(aes(y=get(colnames[4]), colour=colnames[4])) + geom_line(aes(y=get(colnames[5]), colour=colnames[5])) + geom_line(aes(y=get(colnames[6]), colour=colnames[6])) + scale_color_manual(values=color) + ylab("RPM") + labs(colour="") + theme_grey(base_size = 12) + ggtitle(opt$title) + theme(legend.position="right") + ylim(0,1.5)
#dev.off()
