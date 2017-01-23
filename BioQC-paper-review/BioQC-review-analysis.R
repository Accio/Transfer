## Script to benchmark the signatures of BioQC using independent data sources and to test their robustness against potential batch effects in the original data
## I thank two anonymous reviewers of the BioQC paper for suggesting the analysis

## Jitao David Zhang <jitao_david.zhang@roche.com>, 23.01.2017

library(ribiosIO)
library(ribiosUtils)
library(ribiosUDIS)
library(ribiosExpression)
library(ribiosMath)
library(BioQC)
library(limma)
library(sva)
library(ribiosMath)
library(ribiosPlot)
library(ribiosAnnotation)
library(affy)
set.seed(1887)

infile <- function(x) file.path("data", x)
figfile <- function(x) file.path("figures", x)

cacheFile <- infile("esetsCache.RData")
if(!loadFile(cacheFile)) {
    ## the GNF dataset (studyId=481, datasetId=3)
    gnf <- getUDISexpression(id="3", idType="datasetId")
    ## the NB dataset (datasetId=104)
    nb <- getUDISexpression(id="104", idType="datasetId")
    ## mouse dataset (datasetId=1214)
    mouseGNF <- getUDISexpression(id="1214", idType="datasetId")
    ## RNAseq atlas (studyId = 2454, datasetId=461), which is actually combined RNASeqAtlas and BodyMap
    atlas <- getUDISexpression(id="361", idType="datasetId")
    save(gnf, nb, mouseGNF, atlas, file=cacheFile)
}

## load BioQC signatures
bioqcSig <- read_gmt_list(system.file("extdata", "exp.tissuemark.affy.roche.symbols.gmt", package="BioQC"))
logOR <- function(logList) {
    lens <- sapply(logList, length)
    stopifnot(ulen(lens)==1)
    res <- rep(FALSE, lens[[1]])
    for(i in seq(along=logList))
        res <- res | logList[[i]]
    return(res)
}
annotateBioQC <- function(gmt) {
    names <- names(gmt)
    glen <- sapply(gmt, function(x) length(x$genes))
    format <- "^(.*)(NR|NGS_RNASEQATLAS)_([0-9\\.]*)_([0-9\\.]*)$"
    origTissue <- gsub("_$", "", gsub(format, "\\1", names))
    tissue <- trim(gsub("_", " ", origTissue))
    type <- gsub(format, "\\2", names)
    gini <- as.numeric(gsub(format, "\\3", names))
    rank <- as.numeric(gsub(format, "\\4", names))
    res <- data.frame(Name=names,
                      OrigTissue=origTissue,
                      Tissue=tissue,
                      DataSource=type,
                      Gini=gini,
                      Rank=rank,
                      GeneCount=glen)
    return(res)
}
bioqcSigAnno <- annotateBioQC(bioqcSig)

fixConc <- function(x) {
    lapply(x, function(xx) {
               if(is.null(xx)) return(NULL)
               unique(unlist(strsplit(xx, "///")))
           })
}

histMedian <- function(eset, ...) {
    x <- apply(exprs(eset), 2, median)
    hist(x, ...)
}

medianNorm <- function(eset) {
    mat <- exprs(eset)
    colMedian <- apply(mat, 2, median)
    ##meanMed <- mean(colMedian)
    multiSignal <- 50
    normMat <- matrix(rep(colMedian, nrow(mat)), nrow=nrow(mat), byrow=TRUE)
    norm <- mat/normMat*multiSignal
    exprs(eset) <- norm
    return(eset)
}


eset2gini <- function(eset, method=c("Gini", "Shannon"), specThr=0.7, rankThr=3, expThr=1) {
    method <- match.arg(method)
    eset <- medianNorm(eset)
    tissueType <- eset$Tissue.or.Cell.type
    matLists <- tapply(1:ncol(eset), tissueType, function(x) {
                      if(length(x)==1)
                          return(exprs(eset)[,x])
                      return(apply(exprs(eset)[,x], 1, median))
                  })
    mat <- do.call(cbind, matLists)
    if(method=="Gini") {
        matSpec <- gini(mat)
    } else {
        matSpec <- apply(mat, 1, entropy)
    }
    nc <- ncol(mat)
    matRank <- t(apply(mat, 1, function(x) nc-rank(x)+1))
    genes <- fData(eset)$Gene.Symbol
    sigGenes <- lapply(1:ncol(matRank), function(i) {
                           isHit <- matRank[,i]<=rankThr & matSpec>specThr & mat[, i]>expThr
                           hitGenes <- setdiff(as.character(genes[isHit]), c("-", "", " "))
                           return(hitGenes)
                       })
    names(sigGenes) <- levels(tissueType)
    sigGenes <- fixConc(sigGenes)
    return(sigGenes)
}

svaEset <- function(eset, mod=NULL, isLinear=TRUE) {
    if(is.null(mod)) {
        mod <- model.matrix(~eset$Tissue.or.Cell.type)
        colnames(mod) <- levels(eset$Tissue.or.Cell.type)
    }
    mat <- exprs(eset)
    if(isLinear) {
        mat <- log2(mat+0.01)
    }
    eset.sva <- sva(mat, mod=mod)
    sv <- eset.sva$sv
    mat.woBatch <- removeBatchEffect(mat, covariates=sv[,1:5], design=mod)
    if(isLinear) {
        mat.woBatch <- 2^mat.woBatch-0.01
    }
    exprs(eset) <- mat.woBatch
    return(eset)
}

annotateMouseGeneSymbols <- function(gs)  {
    humanGs <- annotateGeneSymbols(gs, organism="mouse", ortholog=TRUE)$GeneSymbol
    humanGs <- unique(humanGs)
    res <- humanGs[!is.na(humanGs)]
    return(res)
}
## SVA does not apply to atlas because no biological replicates are available
## atlasGini <- eset2gini(atlas)

jaccard <- function(x,y) length(intersect(x,y))/length(union(x,y))

svaCacheFile <- infile("svaCache.RData")
if(!loadFile(svaCacheFile)) {
    ## gnf
    gnfGini <- eset2gini(gnf, expThr=50)
    gnfSVA <- svaEset(gnf, isLinear=TRUE)
    gnfSvaGini <- eset2gini(gnfSVA, expThr=50)

    ## svaEffect
    gnfSvaJac <- sapply(seq(along=gnfGini), function(i)
        sapply(seq(along=gnfSvaGini), function(j)
            jaccard(gnfGini[[i]], gnfSvaGini[[j]])))
    colnames(gnfSvaJac) <- rownames(gnfSvaJac) <- names(gnfGini)
    biosHeatmap(gnfSvaJac,
                main="GNF dataset", cexRow=0.85, cexCol=0.85, lwid=c(1,5), lhei=c(1,10),
                xlab="Tissue signatures without SVA correction",
                ylab="Tissue signatures with SVA correction",
                Rowv=FALSE, Colv=FALSE, col="blackyellow", color.key.title="Jaccard index")
    ipdf(figfile("GNF-SVA-jaccardHeatmap.pdf"), width=9L, height=9L)
    
    ## it seems that by applying SVA, many more genes are identified
    compactPar()
    plot(sapply(gnfGini, length), sapply(gnfSvaGini, length),
         main="GNF", 
         xlab="Size of signatures [no SVA]",
         ylab="Size of signatures [with SVA]")
    abline(0, 1, lty=2)
    ipdf(figfile("GNF-SVA-scatterPlot.pdf"), width=5L, height=5L)
    
    ## on average ~60% signatures without adjustment are found also in signatures with SVA adjustments
    summary(gnfSvaEff <- sapply(seq(along=gnfGini), function(i) mean(gnfGini[[i]] %in% gnfSvaGini[[i]])))

    ## gnf mouse
    mouseGNFGini <- eset2gini(mouseGNF, expThr=50)
    mouseGNFSVA <- svaEset(mouseGNF, isLinear=TRUE)
    mouseGNFSvaGini <- eset2gini(mouseGNFSVA, expThr=50)
    mouseGNFHumanGini <- lapply(mouseGNFGini, annotateMouseGeneSymbols)
    mouseGNFHumanSvaGini <- lapply(mouseGNFSvaGini, annotateMouseGeneSymbols)

    mouseGNFSvaJac <- sapply(seq(along=mouseGNFGini), function(i)
        sapply(seq(along=mouseGNFSvaGini), function(j)
            jaccard(mouseGNFGini[[i]], mouseGNFSvaGini[[j]])))
    colnames(mouseGNFSvaJac) <- rownames(mouseGNFSvaJac) <- names(mouseGNFGini)
    biosHeatmap(mouseGNFSvaJac,
                main="MOUSEGNF dataset", cexRow=0.85, cexCol=0.85, lwid=c(1,5), lhei=c(1,10),
                xlab="Tissue signatures without SVA correction",
                ylab="Tissue signatures with SVA correction",
                Rowv=FALSE, Colv=FALSE, col="blackyellow", color.key.title="Jaccard index")
    ipdf(figfile("MOUSEGNF-SVA-jaccardHeatmap.pdf"), width=9L, height=9L)
    
    
    ## NB
    nbGini <- eset2gini(nb, expThr=50)
    nbSVA <- svaEset(nb, isLinear=TRUE)
    nbSvaGini <- eset2gini(nbSVA, expThr=50)

    nbSvaJac <- sapply(seq(along=nbGini), function(i)
        sapply(seq(along=nbSvaGini), function(j)
            jaccard(nbGini[[i]], nbSvaGini[[j]])))
    colnames(nbSvaJac) <- rownames(nbSvaJac) <- names(nbGini)
    biosHeatmap(nbSvaJac,
                main="NB dataset", cexRow=0.85, cexCol=0.85, lwid=c(1,5), lhei=c(1,10),
                xlab="Tissue signatures without SVA correction",
                ylab="Tissue signatures with SVA correction",
                Rowv=FALSE, Colv=FALSE, col="blackyellow", color.key.title="Jaccard index")
    ipdf(figfile("NB-SVA-jaccardIndexHeatmap.pdf"), width=9L, height=9L)
    
    ## it seems that by applying SVA, less genes are identified in NB
    plot(sapply(nbGini, length), sapply(nbSvaGini, length),
         main="NB", 
         xlab="Size of signatures [no SVA]",
         ylab="Size of signatures [with SVA]")
    abline(0, 1, lty=2)
    ipdf(figfile("NB-SVA-scatterPlot.pdf"), width=5L, height=5L)
    
    ## on average ~60% signatures are captured
    summary(nbSvaEff <- sapply(seq(along=nbGini), function(i) mean(nbGini[[i]] %in% nbSvaGini[[i]])))

    gnaSvaCoef <- sapply(1:nrow(gnf), function(i) cor(log10(exprs(gnf)[i,]), log10(exprs(gnfSVA)[i,])))
    nbSvaCoef <- sapply(1:nrow(nb), function(i) cor(log10(exprs(nb)[i,]), log10(exprs(nbSVA)[i,])))
    
    save(gnfSVA, gnfGini, gnfSvaGini, gnfSvaJac,
         mouseGNFSVA, mouseGNFGini, mouseGNFSvaGini, mouseGNFSvaJac,mouseGNFHumanGini, mouseGNFHumanSvaGini,
         nbSVA,nbGini, nbSvaGini, nbSvaJac,
         file=svaCacheFile)
}

## merge into non-redundant genesets
fixNames <- function(names) {
    names <- as.character(names)
    names <- gsub(",", "", names)
    names <- gsub(" +", " ", names)
    return(names)
}
fixSignatureNames <- function(gs) {
    gsNames <- names(gs)
    fixedNames <- fixNames(gsNames)
    matchedNames <- bioqcTissues[match(toupper(fixedNames), toupper(bioqcTissues))]
    data.frame(Name=fixedNames,
               FixedNames=matchedNames)
}
bioqcTissues <- as.character(subset(bioqcSigAnno, DataSource=="NR")$Tissue)
nbGiniFix <- fixSignatureNames(nbGini)
gnfGiniFix <- fixSignatureNames(gnfGini)

##writeLines(bioqcTissues, infile("BioQC-tissueList.txt"))
##writeMatrix(nbGiniFix, infile("NB-tissueFix.txt"), row.names=FALSE)
##writeMatrix(gnfGiniFix, infile("GNF-tissueFix.txt"), row.names=FALSE)
makeLabel <- function(df) {
    label <- paste(gsub(" +", "_", trim(df$FixedNames)), "_NR_0.7_3", sep="")
    df$Label <- label
    return(df)
}
nbGiniAnno <- makeLabel(readMatrix("data/NB-tissueFixed.txt", row.names=FALSE, as.matrix=FALSE))
gnfGiniAnno <- makeLabel(readMatrix("data/GNF-tissueFixed.txt", row.names=FALSE, as.matrix=FALSE))
ulen(uniqueFixedTerms <- setdiff(c(as.character(nbGiniAnno[,2]), as.character(gnfGiniAnno[,2])), ""))
mergeSignatures <- function(gsParent, gsParentAnno,
                            gs1, gs1Anno,
                            gs2, gs2Anno,
                            mode=c("replace", "append")) {

    mode <- match.arg(mode)
    ## combine
    gsAnno <- rbind(gs1Anno, gs2Anno)
    
    ## put RNA seq data aside
    isRNAseq <- with(gsParentAnno, DataSource!="NR")
  
    ## fixed names
    allFixedNames <- setdiff(gsAnno$FixedNames,"")

    ## put NR-unique genesets aside
    isNRunique <- !isRNAseq & !gsParentAnno$Tissue %in% allFixedNames

    ## to be replaced
    isToReplace <- !isRNAseq & !isNRunique
    
    ## merge non-redandant genesets
    fixedGeneSets <- lapply(allFixedNames, function(name) {
                                if(mode=="replace") {
                                    res <- c()
                                } else {
                                    ind <- match(name, gsParentAnno$Tissue)
                                    if(is.na(ind)) {
                                        message("New tissue:", name)
                                        res <- c()
                                    } else {
                                        res <- gsParent[[ind]]
                                    }
                                }
                                if(name %in% gs1Anno$FixedNames) {
                                    ind <- mmatch(name, gs1Anno$FixedNames)[[1]]
                                    genes <- unique(unlist(gs1[ind]))
                                    res <- c(res, genes)
                                }
                                if(name %in% gs2Anno$FixedNames) {
                                    ind <- mmatch(name, gs2Anno$FixedNames)[[1]]
                                    genes <- unique(unlist(gs2[ind]))
                                    res <- c(res, genes)
                                }
                                res <- setdiff(unique(res), "-")
                                return(res)
                            })
    allFixedLabels <- matchColumn(allFixedNames, gsAnno, "FixedNames")$Label
    names(fixedGeneSets) <- allFixedLabels
    fixedGeneSetsAnno <- data.frame(Name=allFixedLabels,
                                    OrigTissue=allFixedNames,
                                    Tissue=allFixedNames,
                                    DataSource="microarray",
                                    Gini=0.7, Rank=3, GeneCount=sapply(fixedGeneSets, length))

    resSignatures <- c(gsParent[isRNAseq], gsParent[isNRunique], fixedGeneSets)
    resAnno <- rbind(gsParentAnno[isRNAseq,], gsParentAnno[isNRunique,], fixedGeneSetsAnno)
    res <- list(signatures=resSignatures,
                annotation=resAnno)
    return(res)
}

bioqcSigGenes <- lapply(bioqcSig, function(x) x$genes)

## without SVA
giniSigAnno <- mergeSignatures(bioqcSigGenes, bioqcSigAnno,
                               nbGini,nbGiniAnno,
                               gnfGini, gnfGiniAnno, mode="replace")

giniSigGenesets <- giniSigAnno$signatures
giniSigAnno  <- giniSigAnno$annotation

## with SVA
giniSvaSigAnno <- mergeSignatures(bioqcSigGenes, bioqcSigAnno,
                                  nbSvaGini,nbGiniAnno,
                                  gnfSvaGini, gnfGiniAnno, mode="replace")

giniSvaSigGenesets <- giniSvaSigAnno$signatures
giniSvaSigAnno  <- giniSvaSigAnno$annotation

## with SVA and appending
giniSvaAppendSigAnno <- mergeSignatures(bioqcSigGenes, bioqcSigAnno,
                                        nbSvaGini,nbGiniAnno,
                                        gnfSvaGini, gnfGiniAnno, mode="append")
giniSvaAppendSigGenesets <- giniSvaAppendSigAnno$signatures
giniSvaAppendSigAnno  <- giniSvaAppendSigAnno$annotation

## compare gene list generated by me and laura - the difference is likely caused by GNF mouse data
meLaura <- (merge(giniSigAnno, bioqcSigAnno, "Name", suffix=c(".me", ".laura")))
with(meLaura, plot(GeneCount.laura~GeneCount.me, log="xy", pch=21, bg=ifelse(DataSource.me=="microarray", "black", "lightgray"),
                   xlab="R-implementation",
                   ylab="Laura's original implementation",
                   main="Size of tissue signatures by implementation"))
ipdf(figfile("tissue-signature-size-byImplementation.pdf"))

## compare with and without SVA
worwoSVA <- merge(giniSigAnno, giniSvaSigAnno, "Name", suffix=c(".noSVA", ".SVA"))
with(worwoSVA, plot(GeneCount.noSVA~GeneCount.SVA, log="xy", pch=21, bg=ifelse(DataSource.noSVA=="microarray", "black", "lightgray"),
                    xlab="R-implementation with SVA",
                    ylab="R-implementation without SVA",
                    main="Size of tissue signatures by SVA"))
ipdf(figfile("tissue-signature-size-bySVA.pdf"))

## compare original and appended
oriVsSVAappend <- merge(bioqcSigAnno, giniSvaAppendSigAnno, "Name", suffix=c(".origin", ".SVAappend"))
with(oriVsSVAappend, plot(GeneCount.SVAappend~GeneCount.origin, log="xy", pch=21,
                          bg=ifelse(DataSource.SVAappend=="microarray", "black", "lightgray"),
                    ylab="R-implementation with SVA hits appended",
                    xlab="Original implementation",
                    main="Size of tissue signatures by SVA"))
ipdf(figfile("tissue-signature-SVAappend.pdf"))

## output
write_gmt(giniSigGenesets, infile("GiniTissueSignatures-Rimp.gmt"), description="Tissue signatures derived from the R implementation of Gini index")
write_gmt(giniSvaSigGenesets,
          infile("GiniTissueSignatures-Rimp-SVA.gmt"), description="Tissue signatures derived from the R implementation of Gini index, with SVA correction")
write_gmt(bioqcSigGenes,
          infile("GiniTissueSignatures-BioQC.gmt"))
write_gmt(giniSvaAppendSigGenesets,
          infile("GiniTissueSignatures-BioQC-SVAappend.gmt"))

## only output replaced ones (107)
repNames <- mintersect(rownames(subset(giniSigAnno, DataSource=="microarray")),
                   names(giniSigGenesets), names(bioqcSigGenes))
giniSigRep <- giniSigGenesets[repNames]
giniSvaSigRep <- giniSvaSigGenesets[repNames]
bioqcSigRep <- bioqcSigGenes[repNames]
bioqcSvaSigRep <- giniSvaAppendSigGenesets[repNames]

write_gmt(giniSigRep, 
          infile("RepSubset-GiniTissueSignatures-Rimp.gmt"), description="Tissue signatures derived from the R implementation of Gini index")
write_gmt(giniSvaSigRep,
          infile("RepSubset-GiniTissueSignatures-Rimp-SVA.gmt"), description="Tissue signatures derived from the R implementation of Gini index, with SVA correction")
write_gmt(bioqcSigRep,
          infile("RepSubset-GiniTissueSignatures-BioQC.gmt"), description="BioQC signatures")
write_gmt(bioqcSvaSigRep,
          infile("RepSubset-GiniTissueSignatures-BioQC-SVAappend.gmt"), description="BioQC signatures with addition of genes identified by the SVA procedure")

## test different variants of replaced signatures using GTEx data
library(readr)
##readGTEx <- function() {
##    gct <- read_gct_matrix("data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct")
##    annotation <- readMatrix("data/GTEx-UDISDataSetID5681-sampleAnnotation.txt", as.matrix=FALSE, row.names=FALSE)
##    annotation.tissue <- readMatrix("data/GTEx-v6_UDIS.txt", as.matrix=FALSE, row.names=FALSE)
##    ids <- colnames(gct)
##    annotation.matched <- matchColumn(ids, annotation, "Experiment name")
##    annotation.tissue.matched <- matchColumn(ids, annotation.tissue, "SAMPID")
##    
##    colData <- cbind(annotation.matched, annotation.tissue.matched)
##    rownames(colData) <- colData[["Experiment name"]]
##    
##    ## filter rows with NA
##    ## RIN: RNA quality index
##    ## SMTS: tissue
##    validInds <- which(!apply(is.na(colData[,c("RIN", "Gender", "SMTS")]), 1, sum) >= 1)
##    gct <- gct[,validInds]
##    colData <- colData[validInds,]
##    
##    rowData <- annotateEnsembl(gsub("\\.[0-9]*$", "", rownames(gct)))
##    rownames(rowData) <- rownames(gct)
##    eset <- new("ExpressionSet",
##                exprs=gct,
##                featureData=new("AnnotatedDataFrame", rowData),
##                phenoData=new("AnnotatedDataFrame", colData))
##    return(eset)
##}
##
##gtexFile <- "data/GTEx-ExpressionSet.RData"
##if(!loadFile(gtexFile)) {
##    gtex <- readGTEx()
##    save(gtex, file=gtexFile)
##}

## run BioQC
test <- wmwTest(gtex, 
giniSigRep, 
          infile("RepSubset-GiniTissueSignatures-Rimp.gmt"), description="Tissue signatures derived from the R implementation of Gini index")
write_gmt(giniSvaSigRep,
          infile("RepSubset-GiniTissueSignatures-Rimp-SVA.gmt"), description="Tissue signatures derived from the R implementation of Gini index, with SVA correction")
write_gmt(bioqcSigRep,
          infile("RepSubset-GiniTissueSignatures-BioQC.gmt"), description="BioQC signatures")
write_gmt(bioqcSvaSigRep,
          infile("RepSubset-GiniTissueSignatures-BioQC-SVAappend.gmt"), description="BioQC signatures with addition of genes identified by the SVA procedure")

## test different variants of replaced signatures using GTEx data
library(readr)
readGTEx <- function() {
    gct <- read_gct_matrix("data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct")
    annotation <- readMatrix("data/GTEx-UDISDataSetID5681-sampleAnnotation.txt", as.matrix=FALSE, row.names=FALSE)
    annotation.tissue <- readMatrix("data/GTEx-v6_UDIS.txt", as.matrix=FALSE, row.names=FALSE)
    ids <- colnames(gct)
    annotation.matched <- matchColumn(ids, annotation, "Experiment name")
    annotation.tissue.matched <- matchColumn(ids, annotation.tissue, "SAMPID")
    
    colData <- cbind(annotation.matched, annotation.tissue.matched)
    rownames(colData) <- colData[["Experiment name"]]
    
    ## filter rows with NA
    ## RIN: RNA quality index
    ## SMTS: tissue
    validInds <- which(!apply(is.na(colData[,c("RIN", "Gender", "SMTS")]), 1, sum) >= 1)
    gct <- gct[,validInds]
    colData <- colData[validInds,]
    
    rowData <- annotateEnsembl(gsub("\\.[0-9]*$", "", rownames(gct)))
    rownames(rowData) <- rownames(gct)
    eset <- new("ExpressionSet",
                exprs=gct,
                featureData=new("AnnotatedDataFrame", rowData),
                phenoData=new("AnnotatedDataFrame", colData))
    return(eset)
}

gtexFile <- "data/GTEx-RPKM-ExpressionSet.RData"
if(!loadFile(gtexFile)) {
    library(edgeR)
    humanAnno <- gtiTaxAnnotation(9606)
    humanGenes <- with(humanAnno, GeneID[GeneType=="protein-coding"])
    gtex <- readGTEx()
    gtex <- gtex[!is.na(fData(gtex)$GeneID) & fData(gtex)$GeneID %in% humanGenes,]
    len <- querydb("SELECT * FROM GA_EXONUNIONLENGTH WHERE TAX_ID='9606'", db="bin", user="genome", password="genome")
    ## transform counts to rpkm
    gtexLen <- matchColumn(fData(gtex)$GeneID, len, "GENE_ID")$LEN
    hasNoLen <- is.na(gtexLen)
    gtex <- gtex[!hasNoLen,]; gtexLen <- gtexLen[!hasNoLen]
    gtexCounts <- exprs(gtex)
    gtexRPKM <- rpkm(gtexCounts, gene.length=gtexLen)
    exprs(gtex) <- gtexRPKM

    save(gtex, gtexCounts, file=gtexFile)
}

## run BioQC
bestMatch <- function(matrix, first=3) {
    decRanks <- apply(matrix, 2, function(x) nrow(matrix)-rank(x, ties.method="first")+1)
    isFirsts <- apply(decRanks, 2, function(x) which(x<=first))
    if(first==1)
        isFirsts <- matrix(isFirsts, nrow=1, byrow=TRUE)
    firstScores <- lapply(1:ncol(matrix), function(i) {
                              scores <- matrix[,i][isFirsts[,i]]
                              return(sort(scores, decreasing=TRUE))
                          })
    firstRes <- data.frame(Rank=rep(1:first, ncol(matrix)),
                           Sample=rep(colnames(matrix), each=first),
                           Tissue=as.vector(sapply(firstScores, names)),
                           Score=unlist(firstScores))
    return(firstRes)
}
gtexBestMatch <- function(matrix, first=3) {
    annoTissue <- gtex$`Tissue or Cell type`
    bm <- bestMatch(matrix, first=first)
    bm$AnnoTissue <- rep(annoTissue, each=first)
    return(bm)
}
bestHeatmap <- function(bests) {
    biosHeatmap(with(subset(bests, Rank==1),
                     table(Tissue, AnnoTissue)),
                Rowv=FALSE, Colv=FALSE, cexRow=1.2, cexCol=1.2, col="blackyellow", zlim=c(0,100))
}

gtexTissueAvgScores <- function(scores) {
    tissues <- factor(gtex$`Tissue or Cell type`)
    stopifnot(ncol(scores)==ncol(gtex))
    resList <- tapply(1:ncol(scores), tissues, function(i)
        rowMedians(scores[,i]))
    res <- do.call(cbind, resList)
    colnames(res) <- levels(tissues)
    rownames(res) <- rownames(scores)
    return(res)
}
gtexBioQCFile <- "data/GTEx-BioQC.RData"
if(!loadFile(gtexBioQCFile)) {
    giniSigRepGtex <- wmwTest(gtex, readGmt(infile("RepSubset-GiniTissueSignatures-Rimp.gmt")), valType="abs.log10p.greater")
    giniSvaSigRepGtex <- wmwTest(gtex, readGmt(infile("RepSubset-GiniTissueSignatures-Rimp-SVA.gmt")), valType="abs.log10p.greater")
    bioqcSigRepGtex <- wmwTest(gtex, readGmt(infile("RepSubset-GiniTissueSignatures-BioQC.gmt")), valType="abs.log10p.greater")
    bioqcSvaSigRepGtex <- wmwTest(gtex, readGmt(infile("RepSubset-GiniTissueSignatures-BioQC-SVAappend.gmt")), valType="abs.log10p.greater")

    giniSigRepBests <- gtexBestMatch(giniSigRepGtex, first=5)
    giniSvaSigRepBests <- gtexBestMatch(giniSvaSigRepGtex, first=5)
    bioqcSigRepBests <- gtexBestMatch(bioqcSigRepGtex, first=5)
    bioqcSvaSigRepBests <- gtexBestMatch(bioqcSvaSigRepGtex, first=5)

    giniSigRepAvg <- gtexTissueAvgScores(giniSigRepGtex)
    giniSvaSigRepAvg <- gtexTissueAvgScores(giniSvaSigRepGtex)
    bioqcSigRepAvg <- gtexTissueAvgScores(bioqcSigRepGtex)
    bioqcSvaSigRepAvg <- gtexTissueAvgScores(bioqcSvaSigRepGtex)

    writeMatrix(giniSigRepAvg, infile("GiniTissueSignatures-Rimp-GtexAvgScores.txt"))
    writeMatrix(giniSvaSigRepAvg, infile("GiniTissueSignatures-Rimp-SVA-GtexAvgScores.txt"))
    writeMatrix(bioqcSigRepAvg, infile("GiniTissueSignatures-BioQC-GtexAvgScores.txt"))
    writeMatrix(bioqcSvaSigRepAvg, infile("GiniTissueSignatures-BioQC-SVAappend-GtexAvgScores.txt"))
    
    save(giniSigRepGtex, giniSvaSigRepGtex,
         bioqcSigRepGtex, bioqcSvaSigRepGtex,
         giniSigRepBests, giniSvaSigRepBests,
         bioqcSigRepBests, bioqcSvaSigRepBests,
         giniSigRepAvg, giniSvaSigRepAvg,
         bioqcSigRepAvg, bioqcSvaSigRepAvg,
         file=gtexBioQCFile)
}

bestHeatmap(giniSigRepBests)
bestHeatmap(giniSvaSigRepBests)
bestHeatmap(bioqcSigRepBests)
bestHeatmap(bioqcSvaSigRepBests)

strongHeatmap <- function(mat, min=5, ...) {
    isVal <- apply(mat, 1, function(x) max(x)>=min)
    mat <- mat[isVal,]
    biosHeatmap(mat, ...)
}
biosHeatmap(bioqcSigRepAvg, col="blackyellow", zlim=c(0, 50), main="no SVA")
ipdf(figfile("bioqcSigRepAvg.pdf"))
biosHeatmap(bioqcSvaSigRepAvg, col="blackyellow", zlim=c(0, 50), main="SVA")
ipdf(figfile("bioqcSvaSigRepAvg.pdf"))

writeMatrix(with(subset(giniSigRepBests, Rank==1), table(Tissue, AnnoTissue)), "data/giniSig-conf.txt")
writeMatrix(with(subset(giniSvaSigRepBests, Rank==1), table(Tissue, AnnoTissue)), "data/giniSvaSig-conf.txt")
writeMatrix(with(subset(bioqcSigRepBests, Rank==1), table(Tissue, AnnoTissue)), "data/bioqc-conf.txt")
writeMatrix(with(subset(bioqcSvaSigRepBests, Rank==1), table(Tissue, AnnoTissue)), "data/bioqc-sva-conf.txt")

bioqcSigRepAvgFirsts <- bestMatch(bioqcSigRepAvg, first=3)
bioqcSvaSigRepAvgFirsts <- bestMatch(bioqcSvaSigRepAvg, first=3)
bioqcAvgFirstsMerge <- merge(bioqcSigRepAvgFirsts, bioqcSvaSigRepAvgFirsts,
                             by=c("Rank", "Sample"), suffix=c(".noSVA", ".SVA"))
writeMatrix(bioqcAvgFirstsMerge, "data/bioqc-avgFirst3.txt")


## no sva versus sva
bioqcSvaMerge <- merge(bioqcSigRepBests,
                       bioqcSvaSigRepBests,
                       by=c("Sample", "AnnoTissue", "Rank"),
                       suffix=c(".noSVA", ".SVA"))
bioqcSvaMergeCons <- subset(bioqcSvaMerge, as.character(Tissue.noSVA)==as.character(Tissue.SVA))
library(lattice)
trellis.par.set(compactTrellis())
xyplot(Score.SVA~Score.noSVA, group=AnnoTissue , data=bioqcSvaMergeCons,
       xlab="BioQC score [without SVA]", ylab="BioQC score [with SVA]",
       main="BioQC applied to GTEx expression dataset", xlim=c(-5, 160), ylim=c(-5, 160), 
       scales=list(alternating=1L, tck=c(1,0),
           x=list(at=seq(0,150,50)), y=list(at=seq(0,150,50))), pch=16, 
       abline=list(0, 1))
ipdf("figures/signatureSVA-GTEx-improvement.pdf", width=6L, height=6L)

## further merge with giniSva
giniMerge <- merge(giniSigRepBests,
                   giniSvaSigRepBests,
                   by=c("Sample", "AnnoTissue", "Rank"),
                   suffix=c(".Rgini", ".RginiSVA"))

allMerge <- merge(bioqcSvaMerge, giniMerge, by=c("Sample", "AnnoTissue", "Rank"))
allMergeCons <- subset(allMerge, as.character(Tissue.noSVA)==as.character(Tissue.SVA)
                       & as.character(Tissue.noSVA)==as.character(Tissue.Rgini) &
                           as.character(Tissue.noSVA)==as.character(Tissue.RginiSVA))
pairs(allMergeCons[,c("Score.noSVA", "Score.SVA", "Score.Rgini", "Score.RginiSVA")],
      pch=16)

## output GNF mouse
writeGct(mouseGNF, infile("mouseGNF-signalMatrix.gct"))
writeMatrix(pData(mouseGNF), infile("mouseGNF-phenoData.txt"))
writeMatrix(fData(mouseGNF), infile("mouseGNF-featureData.txt"))
