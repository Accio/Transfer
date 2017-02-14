library(parallel)

## logFCgage: convert a matrix of logFC values (or dgeScores, defined as -log10(p)*sign(logFC)) into pathway activity scores using the GAGE algorithm
logFCgage <- function(mat, gscs) {
    res <- mclapply(1:ncol(mat), function(i) {
                      x <- mat[,i]
                      names(x) <- rownames(mat)
                      tbl <- ribiosNGS:::myGage(x, gscs)
                      tbl$Contrast <- colnames(mat)[i]
                      tbl <- tbl[,c("Category", "Contrast", "GeneSet", "NGenes", "Direction", "PValue", "FDR")]
                      return(tbl)
                  })
    tbl <- do.call(rbind, res)
    return(tbl)
}

## ribiosNGS::myGage is available in https://github.com/Accio/ribios/blob/master/ribiosNGS/R/gse.R, Below is a hard copy as of 14.02.2017
myGage <- function(logFC, gsc) {
    genes <- gsGenes(gsc)
    cate <- gsCategory(gsc)
    gdf <- data.frame(geneset=names(genes), category=cate)
    gage.res <- gage::gage(logFC, gsets=genes, ref=NULL, samp=NULL)
    greater <- as.data.frame(gage.res[[1]][,c("stat.mean", "p.val", "q.val", "set.size")])
    less <- as.data.frame(gage.res[[2]][,c("p.val", "q.val")])
    greater$geneset <- rownames(greater)
    less$geneset <- rownames(less)
    greater <- merge(greater, gdf, by="geneset")
    less <- merge(less, gdf, by="geneset")
    res.raw <- merge(greater, less, by=c("geneset","category"),
                     suffix=c(".greater", ".less"))
    direction <- with(res.raw, ifelse(p.val.less<p.val.greater, "Down", "Up"))
    pVal.pmin <- with(res.raw, ifelse(p.val.less<p.val.greater, p.val.less, p.val.greater))
    pVal <- pVal.pmin * 2; pVal[pVal>1] <- 1

    ## TODO: contributing genes

    res <- data.frame(Category=res.raw$category,
                      GeneSet=res.raw$geneset,
                      NGenes=res.raw$set.size,
                      Direction=direction,
                      PValue=pVal,
                      FDR=rep(NA,length(pVal)))
    res <- subset(res, NGenes>=1 & !is.na(PValue) & !is.nan(PValue))
    res$FDR <- p.adjust(res$PValue, "fdr")
    return(res)
}


logFCgage <- function(edgeResult, gscs) {
    geneSymbols <- fData(edgeResult)$GeneSymbol
    logFCs <- lapply(dgeTables(edgeResult), function(x) {
                         res <- x$logFC
                         names(res) <- x$GeneSymbol
                         return(res)
                     })
    erTables <- lapply(logFCs, function(x) myGage(x, gscs))
    
    erTable <- do.call(rbind, erTables)
    erTable$Contrast <- rep(names(logFCs),sapply(erTables, nrow))
    erTable <- putColsFirst(erTable, c("Category", "Contrast", "GeneSet"))
    rownames(erTable) <- NULL
    
    edgeGse <- as(edgeResult, "EdgeGSE")
    edgeGse@geneSets <- gscs
    edgeGse@method <- "gage"
    edgeGse@enrichTables <- erTable
    return(edgeGse)
}

