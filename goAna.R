#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing gene list(s)"),
	make_option(c("-o", "--outDir"), help="output directory to store pdf image file(s)"),
	make_option(c("-e", "--genome"), help="genome for which to perform the analysis"),
    make_option(c("-d", "--geneIdType"), default="SYMBOL", help="input gene id type (default=%default)"),
    make_option(c("-a", "--annotation"), default="GOTERM_BP_ALL", help="annotation type (default=%default)"),
    make_option(c("-p", "--pValue"), default="0.05", help="p-value (defaut=%default)"),
    make_option(c("-m", "--maxClass"), default=20, help="maximum number of go classes to plot (default=%default)"),
    make_option(c("-n", "--minGene"), default=10, help="minimum genes in a go class (default=%default)"),
    make_option(c("-c", "--compareCluster"), help="input file contains multiple gene lists (one per column)", action="store_true"),
    make_option(c("-f", "--formula"), help="input file contains multiple gene lists (column 1: genes; column 2: class)", action="store_true"),
    make_option(c("-l", "--listAnnotation"), help="just list different GO annotation options", action="store_true"),
    make_option(c("-s", "--sessionFile"), default="go_analysis.Rsession", help="output session file. It will be used, if already exist"),
    make_option(c("-w", "--figWidth"), default="10", help="width of output figure"),
    make_option(c("-t", "--figHeight"), default="20", help="height of output figure"),
    make_option(c("-u", "--allowDuplicates"), default=0, help="Plot top -m classes even if some are common between classes")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$inFile) | is.null(opt$genome) | is.null(opt$outDir)) & is.null(opt$listAnnotation)) {
	cat("\nProgram: goAna.R (R script to plot go enrichment analysis)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(RDAVIDWebService))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))

##It may help to represent the p-values of enrichment as a transformation whereby you represent the, as -10*log(P-value)

#cat(sprintf("A: %s", opt$compareCluster));
#cat(sprintf("B: %s", opt$formula));
#q();
if(!is.null(opt$listAnnotation)) {
    david<-DAVIDWebService$new("pundhir@binf.ku.dk")
    idType("org.Hs.eg.db")
    getIdTypes(david)
    #getAllAnnotationCategoryNames(david)
    ## ?DAVIDWebService (for more details)
} else if(!is.null(opt$compareCluster)) {
    ## read input gene list file
    gene <- read.table(opt$inFile, fill=TRUE, sep="\t")

    opt$sessionFile <- sprintf("%s/go_analysis_compareCluster_%s.Rsession", opt$outDir, opt$annotation)
    if(!file.exists(opt$sessionFile)) {
        ## convert gene symbol to entrex gene id
        gene[gene==""]  <- NA
        #geneList <- lapply(gene, function(x) { list <- queryMany(x[!is.na(x)], scopes=opt$geneIdType, fields="entrezgene", species=opt$genome)$entrezgene; list[which(!is.na(list))]; } )
        if(opt$genome=="mouse") {
            geneList <- lapply(gene, function(x) { list <- bitr(x[!is.na(x)], fromType=opt$geneIdType, toType="ENTREZID", annoDb="org.Mm.eg.db")$ENTREZID; list[which(!is.na(list))]; } )
        } else {
            geneList <- lapply(gene, function(x) { list <- bitr(x[!is.na(x)], fromType=opt$geneIdType, toType="ENTREZID", annoDb="org.Hs.eg.db")$ENTREZID; list[which(!is.na(list))]; } )
        }

        if(opt$geneIdType=="ENSEMBLTRANS") {
            geneList <- unique(geneList[,c(2,3)])
        }

        ## compute go enrichment
        #opt$pValue <- 0.9
        if(opt$annotation=="DAVID"){
            results=compareCluster(geneList, fun="enrichDAVID", idType="ENTREZ_GENE_ID", listType="Gene", annotation="GOTERM_BP_ALL", david.user = "pundhir@binf.ku.dk", species=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
        } else if(opt$annotation=="KEGG_PATHWAY") {
            results=compareCluster(geneList, fun="enrichKEGG", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
        } else if(opt$annotation=="DISEASE_ONTOLOGY") {
            results=compareCluster(geneList, fun="enrichDO", pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
        } else if(opt$annotation=="REACTOME_PATHWAY") {
            results=compareCluster(geneList, fun="enrichPathway", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
        } else {
            results=compareCluster(geneList, fun="enrichGO", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
        }

    } else {
        load(opt$sessionFile)
    }
    #outFile <- sprintf("%s/go_analysis_compareCluster_%s.pdf", opt$outDir, opt$annotation)
    #pdf(outFile)
    #plot(results, includeAll=TRUE, showCategory=NULL)
    #dev.off()

    library(ggplot2)
    data <- summary(results)[,c(1,2,3,4,5,6,7,8,10)]
    if(nrow(data)>=1) {
        outFile <- sprintf("%s/go_analysis_compareCluster_%s.pdf", opt$outDir, opt$annotation)
        data$GeneDensity <- apply(data, 1, function(x) as.numeric(unlist(strsplit(x[4],"/"))[1])/as.numeric(unlist(strsplit(x[4],"/"))[2]))
        data$BgDensity <- apply(data, 1, function(x) as.numeric(unlist(strsplit(x[5],"/"))[1])/as.numeric(unlist(strsplit(x[5],"/"))[2]))
        data$fc <- log2(data$GeneDensity/data$BgDensity)
        if(as.numeric(opt$allowDuplicates)==0) {
            data <- subset(data, !duplicated(ID))
        }
        data <- data[order(data$Cluster), ]
        data_sig <- Reduce(rbind, by(data, data["Cluster"], head, n=as.numeric(opt$maxClass)))
        if(as.numeric(opt$allowDuplicates)==1) {
            data_sig <- data[data$Description %in% data_sig$Description,]
            data_sig$Description <-factor(data_sig$Description, levels=unique(data_sig$Description))
        } else { 
            data_sig$Description <-factor(data_sig$Description, levels=data_sig$Description)
        }
        p <- ggplot(data_sig, aes(x = Cluster, y = Description)) +
            #geom_point(aes(colour=pvalue,size=GeneDensity)) +
            geom_point(aes(colour=pvalue,size=Count)) +
            scale_colour_gradient(low="red", high="white") +
            scale_size(range=c(1,10))
        #ggsave(p, dpi=300, height=as.numeric(opt$maxClass)*1.5, width=length(unique(data_sig$Cluster))*1.5, filename=outFile, useDingbats=FALSE)
        #ggsave(p, dpi=300, height=as.numeric(opt$figHeight), width=as.numeric(opt$figWidth), filename=outFile, useDingbats=FALSE)
        ggsave(p, dpi=300, height=20, width=10, filename=outFile, useDingbats=FALSE)
    }

    outFile <- sprintf("%s/go_analysis_compareCluster_%s.xls", opt$outDir, opt$annotation)
    #data <- as.data.frame(summary(results))
    #data$geneIDf <- apply(data, 1, function(y) paste(unlist(lapply(unlist(strsplit(y[9], "/")), function(x) geneList[which(geneList$ENTREZID==x),]$V1)), collapse="/"))
    #write.table(data, file=outFile, sep="\t", quote=F, col.names=T, row.names=F)
    write.table(as.data.frame(summary(results)), file=outFile, sep="\t", quote=F, col.names=T, row.names=F)

    save.session(opt$sessionFile)
} else if(!is.null(opt$formula)) {
    ## read input gene list file
    geneList <- read.table(opt$inFile)

    opt$sessionFile <- sprintf("%s/go_analysis_compareClusterFormula_%s.Rsession", opt$outDir, opt$annotation)
    if(!file.exists(opt$sessionFile)) {
        ## convert gene symbol to entrex gene id
        if(opt$genome=="mouse") {
            geneList_conv <- bitr(geneList$V1, fromType=opt$geneIdType, toType="ENTREZID", annoDb="org.Mm.eg.db")
        } else {
            geneList_conv <- bitr(geneList$V1, fromType=opt$geneIdType, toType="ENTREZID", annoDb="org.Hs.eg.db")
        }

        geneList_conv <- as.data.frame(geneList_conv)
        geneList <- merge(geneList, geneList_conv, by.x="V1", by.y=opt$geneIdType)
        geneList <- geneList[which(!is.na(geneList$ENTREZID)),]

        if(opt$geneIdType=="ENSEMBLTRANS") {
            geneList <- unique(geneList[,c(2,3)])
        }

        ## compute go enrichment
        if(opt$annotation=="DAVID"){
            results=compareCluster(ENTREZID~V2, data=geneList, fun="enrichDAVID", idType="ENTREZ_GENE_ID", listType="Gene", annotation="GOTERM_BP_ALL", david.user = "pundhir@binf.ku.dk", species=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
        } else if(opt$annotation=="KEGG_PATHWAY"){
            results=compareCluster(ENTREZID~V2, data=geneList, fun="enrichKEGG", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        } else if(opt$annotation=="DISEASE_ONTOLOGY"){
            results=compareCluster(ENTREZID~V2, data=geneList, fun="enrichDO", pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        } else if(opt$annotation=="REACTOME_PATHWAY"){
            results=compareCluster(ENTREZID~V2, data=geneList, fun="enrichPathway", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        } else {
            results=compareCluster(ENTREZID~V2, data=geneList, fun="enrichGO", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        }
    } else {
        load(opt$sessionFile)
        #opt$maxClass=10
        #opt$outDir="go_analysis/pu1_pregm_wt_ko/all"
    }

    data <- summary(results)[,c(1,2,3,4,5,6,7,8,10)]
    if(nrow(data)>1) {
        outFile <- sprintf("%s/go_analysis_compareClusterFormula_%s.pdf", opt$outDir, opt$annotation)
        library(ggplot2)
        data$GeneDensity <- apply(data, 1, function(x) as.numeric(unlist(strsplit(x[4],"/"))[1])/as.numeric(unlist(strsplit(x[4],"/"))[2]))
        data$BgDensity <- apply(data, 1, function(x) as.numeric(unlist(strsplit(x[5],"/"))[1])/as.numeric(unlist(strsplit(x[5],"/"))[2]))
        data$fc <- log2(data$GeneDensity/data$BgDensity)
        if(as.numeric(opt$allowDuplicates)==0) {
            data <- subset(data, !duplicated(ID))
        }
        data <- data[order(data$Cluster), ]
        data_sig <- Reduce(rbind, by(data, data["Cluster"], head, n=as.numeric(opt$maxClass)))
        if(as.numeric(opt$allowDuplicates)==1) {
            data_sig <- data[data$Description %in% data_sig$Description,]
            data_sig$Description <-factor(data_sig$Description, levels=unique(data_sig$Description))
        } else {
            data_sig$Description <-factor(data_sig$Description, levels=data_sig$Description)
        }
        p <- ggplot(data_sig, aes(x = Cluster, y = Description)) +
            geom_point(aes(colour=pvalue,size=Count)) +
            scale_colour_gradient(low="red", high="white") +
            scale_size(range=c(5,10)) +
            theme(text=element_text(size=10), axis.text.x=element_text(angle=90))
        #ggsave(p, dpi=300, height=as.numeric(opt$maxClass)*1.5, width=length(unique(data_sig$Cluster))*2, filename=outFile, useDingbats=FALSE)
        #ggsave(p, dpi=300, height=as.numeric(opt$figHeight), width=as.numeric(opt$figWidth), filename=outFile, useDingbats=FALSE)
        ggsave(p, dpi=300, height=10, width=8, filename=outFile, useDingbats=FALSE)
    }

    outFile <- sprintf("%s/go_analysis_compareClusterFormula_%s.xls", opt$outDir, opt$annotation)
    data <- as.data.frame(summary(results))
    data$geneIDf <- apply(data, 1, function(y) paste(unlist(lapply(unlist(strsplit(y[9], "/")), function(x) geneList[which(geneList$ENTREZID==x),]$V1)), collapse="/"))
    write.table(data, file=outFile, sep="\t", quote=F, col.names=T, row.names=F)

    #opt$sessionFile <- sprintf("%s/go_analysis_compareClusterFormula_%s.Rsession", opt$outDir, opt$annotation)
    save.session(opt$sessionFile)
} else {
    ## read input gene list file
    gene <- read.table(opt$inFile)

    opt$sessionFile <- sprintf("%s/go_analysis_list_%s.Rsession", opt$outDir, opt$annotation)
    if(!file.exists(opt$sessionFile)) {
        ## convert gene symbol to entrez gene id
        if(opt$genome=="mouse") {
            gene_conv <- bitr(gene$V1, fromType=opt$geneIdType, toType="ENTREZID", annoDb="org.Mm.eg.db")
        } else {
            gene_conv <- bitr(gene$V1, fromType=opt$geneIdType, toType="ENTREZID", annoDb="org.Hs.eg.db")
        }
        gene_conv <- as.data.frame(gene_conv)
        gene <- merge(gene, gene_conv, by.x="V1", by.y=opt$geneIdType)
        gene <- gene[which(!is.na(gene$ENTREZID)),]
        #row.names(gene) <- gene$entrezgene

        ## compute go enrichment
        if(opt$annotation=="DAVID"){
            ## annotation options: GOTERM_BP_ALL, GOTERM_MF_ALL, GOTERM_CC_ALL, UP_TISSUE, UCSC_TFBS
            results <- enrichDAVID(gene=gene$ENTREZID, idType="ENTREZ_GENE_ID", listType="Gene", annotation="GOTERM_BP_ALL", david.user="pundhir@binf.ku.dk", species=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
        } else if(opt$annotation=="KEGG_PATHWAY"){
            results=enrichKEGG(gene$ENTREZID, organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene))
            #results=enrichKEGG(gene$ENTREZID, organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), minGSSize=as.numeric(opt$minGene), use_internal_data=TRUE)
        } else if(opt$annotation=="DISEASE_ONTOLOGY"){
            results=enrichDO(gene$ENTREZID, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        } else if(opt$annotation=="REACTOME_PATHWAY"){
            results=enrichPathway(gene$ENTREZID, organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        } else if(opt$annotation=="GOTERM_BP_ALL"){
            results=enrichGO(gene$ENTREZID, ont="BP", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        } else if(opt$annotation=="GOTERM_MF_ALL"){
            results=enrichGO(gene$ENTREZID, ont="MF", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        } else if(opt$annotation=="GOTERM_CC_ALL"){
            results=enrichGO(gene$ENTREZID, ont="CC", organism=opt$genome, pvalueCutoff=as.numeric(opt$pValue), qvalueCutoff=as.numeric(opt$pValue))
        }
    } else {
        load(opt$sessionFile)
    }

    save.session(opt$sessionFile)
    ## compute gene enrichment network for up- and down-regulated genes
    if(is.numeric(gene$V2) & nrow(summary(results))>=1) {
        geneList <- structure(gene$V2, names=gene$ENTREZID)
        pdf(sprintf("%s/go_analysis_list_cnetplot_%s.pdf", opt$outDir, opt$annotation), width=20, height=20)
        cnetplot(results, categorySize="pvalue", foldChange=geneList, showCategory = 3)
        dev.off()
    }

    outFile <- sprintf("%s/go_analysis_list_%s.xls", opt$outDir, opt$annotation)
    data <- as.data.frame(summary(results))
    if(nrow(data)>=1) {
        data$geneIDf <- apply(data, 1, function(y) paste(unlist(lapply(unlist(strsplit(y[8], "/")), function(x) gene[which(gene$ENTREZID==x),]$V1)), collapse="/"))
    }
    write.table(data, file=outFile, sep="\t", quote=F, col.names=T, row.names=F)
    pdf(sprintf("%s/go_analysis_list_%s.pdf", opt$outDir, opt$annotation))
    #dotplot(results, x="count", showCategory=opt$maxClass, colorBy="qvalue")
    barplot(results)
    dev.off()

    save.session(opt$sessionFile)
}
