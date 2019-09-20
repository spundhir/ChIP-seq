#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input BED file created using multiIntersectBed (can be stdin; 5th col: id_1,id_2)"),
	make_option(c("-o", "--outFile"), help="output pdf file"),
	make_option(c("-l", "--list"), help="input file contains list (format: id condition; eg. ENSG00000001617 WT)", action="store_true"),
	make_option(c("-t", "--type"), default="ellipses", help="type of plot; ellipses or ChowRuskey. (default: %default)"),
	make_option(c("-s", "--gs"), default="21000", help="genome size meaning total number of genes as background to compute p-values (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outFile)) {
	cat("\nProgram: multiIntersect2venn.R (R script to plot venn diagram from multiIntersectBed results - max 5 classes)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(Vennerable))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library(GeneOverlap))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(pheatmap))

if(opt$inFile=="stdin") {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$inFile)
}
if(is.null(opt$list)) {
    ## In case of following error
    ## brewer.pal minimal value for n is 3, returning requested palette with 3 different levels (uncomment)
    #vec <- as.vector(unique(data[!grepl("[,_]+", data$V5),]$V5))
    vec <- as.vector(unique(data[!grepl("[,]+", data$V5),]$V5))
    data$id <- sprintf("%s_%d_%d", data$V1, data$V2, data$V3)
} else {
    colnames(data) <- c("id", "V5")
    vec <- as.vector(unique(data[!grepl(",", data$V5),]$V5))
}
#lst <- list()
#k <- 1
#for(i in vec) {
#    if(is.null(opt$list)) {
#        l <- data[grep(i, data$V5),]$id
#    } else {
#        l <- data[grep(sprintf("^%s$", i), data$V5),]$id
#    }
#    lst[[k]] <- l
#    k=k+1
#}

if(is.null(opt$list)) {
    lst <- lapply(vec, function(x) data[grep(x, data$V5),]$id)
} else {
    lst <- lapply(vec, function(x) data[grep(x, data$V5),]$id)
    #lst <- lapply(vec, function(x) data[grep(sprintf("^%s$", x), data$V5),]$id)
}
names(lst) <- vec
col <- brewer.pal(length(vec)+1, "Spectral")
col <- col[1:length(vec)]

#gom.obj <- newGOM(lst, lst, opt$gs)
#mat <- getMatrix(gom.obj, name="pval")

if(length(names(lst)) <= 3 & opt$type=="ellipses") {
    #venn.plot <- venn.diagram(lst, fill=col, NULL)
    #pdf(opt$outFile)
    #grid.draw(venn.plot)
    #pheatmap(mat, display_numbers=T, number_format="%.1e", fontsize=10)
    #dev.off()
    Vstem <- Venn(lst)
    pdf(opt$outFile)
    plot(Vstem)
    #pheatmap(mat, display_numbers=T, number_format="%.1e", fontsize=10)
    dev.off()
} else {
    Vstem <- Venn(lst)
    pdf(opt$outFile)
    plot(Vstem, type=opt$type)
    #pheatmap(mat, display_numbers=T, number_format="%.1e", fontsize=10)
    dev.off()
}

## determine names of overlapping elements
#attr(gplots::venn(lst, show.plot=FALSE), "intersections")

#save.session("overlap_venn.session")
