#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--configFile"), help="input configuration file containing bam file information (can be stdin; FORMAT: tpm test.bam)"),
	make_option(c("-o", "--outFile"), help="output file containing read counts"),
    make_option(c("-n", "--genome"), help="genome (mm9, mm10, hg19, hg38)"),
    make_option(c("-j", "--gtfFile"), help="input file containing genomic coordinate of genes in GTF format (optional)"),
    make_option(c("-t", "--featureType"), default="exon", help="specify the feature type (default: %default)"),
    make_option(c("-r", "--attributeType"), default="gene_name", help="specify the feature type (default: %default)"),
    make_option(c("-O", "--allowMultiassign"), action="store_true", help="allow reads assignment to multiple features"),
    make_option(c("-M", "--countMultimapping"), action="store_true", help="also count multi-mapping reads/fragments"),
    make_option(c("-s", "--strandSpecific"), default=0, help="0: unstranded, 1: stranded, 2: reversely stranded (default: %default)"),
    make_option(c("-T", "--processors"), default=1, help="number of processors (default: %default)"),
    make_option(c("-F", "--fraction"), action="store_true", help="If specified, a fractional count 1/n will be generated for each multi-mapping read")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$configFile) | is.null(opt$outFile) | is.null(opt$genome)) {
	cat("\nProgram: motifDynAna.R (R script to compute read count corresponding to input GTF file)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library("Rsubread"))
suppressPackageStartupMessages(library(session))

if(identical(opt$configFile, "stdin")==T) { 
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$configFile)
}

if(is.null(opt$gtfFile)) {
    counts <- featureCounts(data$V2, annot.inbuilt=opt$genome, isGTFAnnotationFile=T, GTF.featureType=opt$featureType, GTF.attrType=opt$attributeType, allowMultiOverlap=!is.null(opt$allowMultiassign), countMultiMappingReads=!is.null(opt$countMultimapping), strandSpecific=opt$strandSpecific, nthreads=opt$processors, fraction=!is.null(opt$fraction))
} else {
    counts <- featureCounts(data$V2, annot.ext=opt$gtfFile, isGTFAnnotationFile=T, GTF.featureType=opt$featureType, GTF.attrType=opt$attributeType, allowMultiOverlap=!is.null(opt$allowMultiassign), countMultiMappingReads=!is.null(opt$countMultimapping), strandSpecific=opt$strandSpecific, nthreads=opt$processors, fraction=!is.null(opt$fraction))
}

all <- cbind(counts$annotation, round(counts$counts))
write.table(all, opt$outFile, sep="\t", quote = F, col.names = T, row.names = F)

save.session("featureCounts.session")
