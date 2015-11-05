#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--tfExprFile"), help="input file containing TF signal or gene expression values for different NFR dynamic classes"),
    make_option(c("-d", "--samplesName"), help="name of samples for which expression has been measured seperated by comma"),
    make_option(c("-c", "--NFRClassCol"), default="5", help="column containing NFR class description (default: %default)"),
    make_option(c("-m", "--mergeRep"), default="TRUE", help="plot mean expression values of replicates (default: %default)"),
    make_option(c("-n", "--minFreq"), default="50", help="minimum frequency of NFRs in each nfr dynamic class (defaut=%default)"),
    make_option(c("-r", "--orderBy"), default="1", help="order bar plot by variable (1) or class (2)  (defaut=%default)"),
	make_option(c("-o", "--outPdfFile"), help="output pdf image file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$tfExprFile) | is.null(opt$samplesName) | is.null(opt$outPdfFile)) {
	cat("\nProgram: nfrDynAnaWithExpr.R (R script to plot TF signal or gene expression with respect to NFR dynamic classes)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(gridExtra))

samplesName <- as.vector(unlist(strsplit(opt$samplesName, ",")))
uniqueSamplesName <- unique(samplesName)
nSamples <- length(samplesName)

data <- read.table(opt$tfExprFile)
ncol <- ncol(data)
end <- ncol
start <- (ncol-as.numeric(nSamples))+1

## take mean expression for replicates
if(opt$mergeRep) {
    k <- 1
    for(i in seq(start, end, by=2)) {
        data[,(ncol+k)] <- apply(data[,c(i:(i+1))], 1, function(x) mean(x))
        k <- k + 1
    }

    colnames(data)[c((ncol+1):ncol(data))] <- uniqueSamplesName
    data.melt <- melt(data[,c(as.numeric(opt$NFRClassCol),(ncol+1):ncol(data))])
} else {
    colnames(data)[c(start:end)] <- uniqueSamplesName
    data.melt <- melt(data[,c(as.numeric(opt$NFRClassCol),start:end)])
}

colnames(data.melt) <- c("class","variable","value")
#head(data.melt)

## correctly sort barplot from lsk, pregm to gmp
if(!is.na(match("TRUE", grepl("granulocytes", samplesName)))) {
    data.melt$variable <- factor(data.melt$variable, levels=c(as.vector(unique(data.melt$variable[grep("lsk", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("pregm", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("gmp", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("granulocytes", data.melt$variable)]))), ordered=TRUE) 
} else if(!is.na(match("TRUE", grepl("pregm", samplesName)))) {
    data.melt$variable <- factor(data.melt$variable, levels=c(as.vector(unique(data.melt$variable[grep("lsk", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("pregm", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("gmp", data.melt$variable)]))), ordered=TRUE)
} else if(!is.na(match("TRUE", grepl("hsc", samplesName))) & is.na(match("TRUE", grepl("nk", samplesName)))) {
    data.melt$variable <- factor(data.melt$variable, levels=c(as.vector(unique(data.melt$variable[grep("hsc", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("clp", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("gmp", data.melt$variable)])), as.vector(unique(data.melt$variable[grep("erya", data.melt$variable)]))), ordered=TRUE)
}

## compute p-values (http://stackoverflow.com/questions/18335209/how-to-plot-additional-statistics-in-boxplot-for-each-group)
dt <- data.table(data.melt)
pval <- dt[, list(pvalue = paste0("pval = ", sprintf("%.3f", summary(aov(value ~ variable))[[1]][["Pr(>F)"]][1]))), by=class]
pval

## compute frequency of NFR for each NFR dynamics class
#class_summary <- as.data.frame(table(data.melt$class))
class_summary <- as.data.frame(table(data[,as.numeric(opt$NFRClassCol)]))
data.melt$freq <- apply(data.melt, 1, function(x) class_summary[which(class_summary$Var1==x[1]),2])
data.melt$class <- sprintf("%s_%s", data.melt$class, data.melt$freq)
data.melt <- data.melt[which(data.melt$freq >= as.numeric(opt$minFreq)),]
if(grepl("_", data.melt$variable)[1]){
    data.melt$tf <- gsub("_.*", "", data.melt$variable)
} else {
    data.melt$tf <- "NA"
}

#data.melt$class <- factor(data.melt$class, levels=c(
#    as.vector(unique(data.melt$class[grep("^pregm,granulocytes_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^pregm,gmp,granulocytes_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^gmp,granulocytes_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^granulocytes_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^lsk,pregm,gmp,granulocytes_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^gmp_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^pregm,gmp_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^pregm_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^lsk_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^lsk,pregm,gmp_", data.melt$class)])), 
#    as.vector(unique(data.melt$class[grep("^lsk,pregm_", data.melt$class)]))), 
#    ordered=TRUE)

p <- list()
nTf <- length(unique(data.melt$tf))
j <- 1
for(i in seq(1, nTf*2, by=2)) {
    sub_data.melt <- data.melt[which(data.melt$tf==unique(data.melt$tf)[j]),]
    if(opt$orderBy==1) {
        p[[i]] <- ggplot(sub_data.melt, aes(factor(variable), log(value))) + geom_boxplot(aes(fill=factor(class))) + ylab("TPM (log)") +
        #xlab("") + 
        facet_grid(.~class) +
        #scale_fill_manual(values=rep("white", length(unique(data.melt$class)))) + 
        #scale_x_discrete(labels=uniqueSamplesName) +
        theme(text=element_text(size=30), axis.text.x=element_text(angle=90)) +
        stat_summary(fun.y=median, geom="line", aes(group=1), color="red") +
        stat_summary(fun.y=median, geom="point") 
        #geom_line(stat="hline", yintercept="median", aes(group=1), color="red", linetype="dashed")
    } else {
        p[[i]] <- ggplot(sub_data.melt, aes(factor(class), log(value))) + geom_boxplot(aes(fill=factor(variable))) + ylab("TPM (log)") +
        #xlab("") + 
        facet_grid(.~variable) +
        #scale_fill_manual(values=rep("white", length(unique(data.melt$variable)))) + 
        #scale_x_discrete(labels=uniqueSamplesName) +
        theme(text=element_text(size=30), axis.text.x=element_text(angle=90)) +
        stat_summary(fun.y=median, geom="line", aes(group=1), color="red") +
        stat_summary(fun.y=median, geom="point") 
        #geom_line(stat="hline", yintercept="median", aes(group=1), color="red", linetype="dashed")
    }

    ## uncomment to order bar plot by gene expression
    #p[[i+1]] <- ggplot(sub_data.melt, aes(x=class, log(value))) +
    p[[i+1]] <- ggplot(sub_data.melt, aes(x=reorder(class, -value, FUN=median), log(value))) +
    geom_boxplot() +
    theme(text=element_text(size=30), axis.text.x=element_text(angle=90)) +
    geom_hline(aes(yintercept=median(log(sub_data.melt$value))), colour="red", linetype="dashed")
    j <- j+1
}
#ggsave(p, dpi=300, height=10, width=40, filename=opt$outPdfFile, useDingbats=FALSE, limitsize=FALSE)
#ggsave(p[[1]], dpi=300, height=10, width=40, filename=opt$outPdfFile, useDingbats=FALSE, limitsize=FALSE)
ggsave(opt$outPdfFile, dpi=300, height=(length(p)*10), width=60, do.call(arrangeGrob, p), limitsize=FALSE)
save.session("test.session")
#warnings()
