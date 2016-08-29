#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inDirs"), help="input directories containing motifStatAna results separated by a comma"),
	make_option(c("-o", "--outDir"), help="output directory to keep pdf image file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inDirs) | is.null(opt$outDir)) {
	cat("\nProgram: motifStatAna.R (R script to plot motif statistics)\n")
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
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

## create output directory
dir.create(opt$outDir)

## parse input directory information
inDirs <- as.vector(unlist(strsplit(opt$inDirs, ",")))
inDirsDes <- gsub(".*/", "", gsub("/$","", inDirs))

## plot histograms of motif location with respect to peaks
all <- NULL
lst <- vector("list", length(inDirs))
for(i in 1:length(inDirs)) {
    df <- read.table(sprintf("%s/REGIONS_INTEREST.hist", inDirs[i]), sep="\t", header = T)
    colnames(df) <- gsub("\\..*", "", colnames(df))
    df$class <- inDirsDes[i]
    lst[[i]] <- df
}
all <- rbind(all, do.call(rbind, lst))
all$Distance <- factor(all$Distance)
pdf(sprintf("%s/motif_hist.pdf", opt$outDir))
for(i in 2:(ncol(all)-1)) {
    all.melt <- melt(all[,c(1,i,ncol(all))])
    print(ggplot(all.melt, aes(Distance, value, fill=factor(variable))) +
    geom_bar(stat="identity") + facet_grid(.~class) +
    ylab(sprintf("Sites per bp per peak (%s)", unique(all.melt$variable))) +
    theme_bw() + theme(legend.position="none") +
    theme(text=element_text(size=15), axis.text.x=element_text(angle = 0, vjust=1)))
}
dev.off()

## plot match score to PWM of motifs  
all <- NULL
all_freq <- NULL
lst <- vector("list", length(inDirs))
lst_freq <- vector("list", length(inDirs))
total <- NULL
for(i in 1:length(inDirs)) {
    df <- read.table(sprintf("%s/REGIONS_INTEREST.find", inDirs[i]), sep="\t", header = T)
    colnames(df) <- gsub("\\..*", "", colnames(df))
    df$class <- inDirsDes[i]
    lst[[i]] <- df

    df <- as.data.frame(t(table(df[,c(4,1)])))
    df$class <- inDirsDes[i]
    lst_freq[[i]] <- df

    total[i] <- countLines(sprintf("%s/REGIONS_INTEREST.bed", inDirs[i]))[1]
}
all <- rbind(all, do.call(rbind, lst))
all_freq <- rbind(all_freq, do.call(rbind, lst_freq))

all.melt <- melt(all[,c(4,6,7)])
p1 <- ggplot(all.melt, aes(class, value)) + geom_boxplot(aes(col=class)) + theme_bw() +
        ylab("Match score to PWM") + theme(legend.position="none") +
        facet_grid(.~Motif) +
        theme(text=element_text(size=15), axis.text.x=element_text(angle = 0, vjust=1))

## plot fraction of peaks that contain the motifs  
df <- unique(all[,c(1,4,7)])
df <- as.data.frame(table(df[,c(3,2)]))
df$total <- rep(total, nrow(df)/length(total))
df$percentage <- (df$Freq/df$total)*100
p2 <- ggplot(df, aes(class, percentage, fill=factor(Motif))) +
        geom_bar(stat="identity") +
        geom_text(aes(label = sprintf("%0.0f", percentage), y = percentage), size = 3) +
        theme_bw() + scale_fill_brewer(palette = "Set3", type = "qualitative") +
        #theme(legend.position="none") +
        theme(text=element_text(size=15), axis.text.x=element_text(angle = 0, vjust=1))

## plot frequency of motifs at each peak 
all.melt <- melt(all_freq[,c(2,3,4)])
p3 <- ggplot(all.melt, aes(class, value)) + geom_boxplot(aes(col=class)) + theme_bw() +
        ylab("Frequency of motifs") +
        facet_grid(.~Motif) + theme(legend.position="none") +
        theme(text=element_text(size=15), axis.text.x=element_text(angle = 0, vjust=1))

g <- arrangeGrob(p1, p2, p3, ncol = 1, nrow = 3)
ggsave(g, dpi=300, height=25,width=25, filename=sprintf("%s/motif_find.pdf", opt$outDir), useDingbats=FALSE)
