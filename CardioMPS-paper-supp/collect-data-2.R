## data mining for the Cardiomyopathy paper

library(GEOquery)
library(ribiosUtils)
library(ribiosIO)
library(limma)
library(ribiosAnnotation)

BASEDIR <- "./data-mining-2-destdir"

gseIDs <- c("GSE82188","GSE54893","GSE68857",
            "GSE71613","GSE71912","GSE63759",
            "GSE67492","GSE65446","GSE63847",
            "GSE54681","GSE82290","GSE52601")
## get all datasets
studyFile <- file.path(BASEDIR, "dataMining2.RData")
if(!loadFile(studyFile)) {
    geoDataSets <- sapply(gseIDs, getGEO, destdir=BASEDIR)
    
    ## NGS studies have not signals
    ## microarray
    isNGS <- sapply(geoDataSets, function(x) nrow(x)==0)
    arrayDatasets <- geoDataSets[!isNGS]
    (arrayIDs <- gseIDs[!isNGS])
##    [1] "GSE82188" "GSE54893" "GSE68857" "GSE63759" "GSE67492" "GSE63847" "GSE54681"
##    [8] "GSE82290" "GSE52601"
    names(arrayDatasets) <- gseIDs[!isNGS]
    
    ## NGS datasets
    ngsIDs <- c("GSE71613", "GSE71912", "GSE65446", "GSE64391")
    ngsDs <- sapply(ngsIDs, getGEO, destdir=BASEDIR)
    names(ngsDs) <- ngsIDs
    lapply(ngsIDs, function(x) getGEOSuppFiles(x,
                                           makeDirectory=TRUE,
                                           baseDir=BASEDIR))

    ## decompress using 'find'

    save(arrayDatasets, ngsGSEs,
         ngsIDs, ngsDs,
         file=studyFile)
}

outfile <- function(x) file.path(file.path(BASEDIR, 'data-for-DGE'), x)
varpd <- function(eset) {
    df <- removeInvarCol(pData(eset))
    toRemove <- c("geo_accession",
                  grep("supplementary_file", colnames(df), value=TRUE),
                  grep("relation", colnames(df), value=TRUE)
                  )
    df[, !colnames(df) %in% toRemove]
}
outfileid <- function(fmt, id) outfile(sprintf(fmt, id))
exportStudy <- function(eset, id, log2=TRUE, design, contrast, featureID=NULL) {
    fd <- fData(eset)
    pd <- pData(eset)
    geneExp <- exprs(eset)
    if(log2)
        geneExp <- log2(geneExp)
    write_gct(fd, file=outfileid("%s-originalFeatureAnnotation.txt", id))
    if(is.null(featureID))
        featureID <- rownames(geneExp)
    stopifnot(length(featureID) == nrow(eset))
    features <- annotateAnyIDs(featureID, orthologue=TRUE, multiOrth=FALSE)

    features$Input <- rownames(geneExp)
    writeMatrix(features, file=outfileid("%s-featureAnnotation.txt", id), row.names=FALSE)
    writeMatrix(pd, file=outfileid("%s-sampleAnnotation.txt", id))
    write_gct(geneExp, file=outfileid("%s-expression.gct", id))
    writeMatrix(design, file=outfileid("%s-designMatrix.txt", id))
    writeMatrix(contrast, file=outfileid("%s-contrastMatrix.txt", id))
}
exportArray <- function(id, log2=TRUE, design, contrast, featureID=NULL) {
    stopifnot(id %in% names(arrayDatasets))
    eset <- arrayDatasets[[id]]
    exportStudy(eset, id, log2=log2, design=design, contrast=contrast, featureID=featureID)
}

##----------------------------------------##
## Old studies
##----------------------------------------##
loadFile("data/cardiomyopathy-meta.RData")

exportStudy(wittchen, id="GSE4172", log2=FALSE, design=wittchen.model, contrast=wittchen.contrast[,"DCMi",drop=FALSE])
exportStudy(gaertner, id="GSE29819", log2=FALSE, design=gaertner.model, contrast=gaertner.contrast[,c(1,2,4,5)])
exportStudy(kittleson, id="GSE1869", log2=FALSE, design=kittleson.model, contrast=kittleson.contrast[,"ICM", drop=FALSE])
exportStudy(kuner, id="GSE3586", log2=FALSE, design=kuner.model, contrast=kuner.contrast,
            featureID=as.character(fData(kuner)$GeneID))
exportStudy(barth, id="GSE3585", log2=FALSE, design=barth.model, contrast=barth.contrast)
exportStudy(harding, id="GSE16909", log2=FALSE, design=harding.model, contrast=harding.contrast,
            featureID=as.character(fData(harding)$Entrez_Gene_ID))
exportStudy(glynjones, id="GSE5606", log2=FALSE, design=glynjones.model, contrast=glynjones.contrast)
## exportStudy(li, id="GSE51483", log2=FALSE, design=li.model, contrast=li.contrast)

##----------------------------------------##
## NEW STUDIES
##----------------------------------------##

## GSE81288: calrecticulin-induced DCM
gse82188 <- arrayDatasets[["GSE82188"]]
gse82188.design <- with(pData(gse82188), model.matrix(~characteristics_ch1.1))
colnames(gse82188.design) <- c("Control","CalreticulinDCM")
gse82188.contrast <- makeContrasts(contrasts="CalreticulinDCM", levels=gse82188.design)
exportArray("GSE82188", log2=TRUE, design=gse82188.design, contrast=gse82188.contrast)

## GSE54893
g893 <- arrayDatasets[["GSE54893"]]
head(exprs(g893)) ## log=TRUE
g893.pd <- varpd(g893)
g893.pd$MCAT <- with(g893.pd, grepl("MCAT", characteristics_ch1.1))
g893.pd$AZT <- with(g893.pd, grepl("AZT", characteristics_ch1.3))
g893.design <- with(g893.pd, model.matrix(~MCAT*AZT))
colnames(g893.design) <- c("Baseline", "MCAT", "AZT", "MCAT_AZT")
g893.contrast <- makeContrasts("MCATwithAZT_vs_MCAT"="AZT+MCAT_AZT",
                              levels=colnames(g893.design))
exportArray("GSE54893", log2=TRUE, design=g893.design, contrast=g893.contrast)

## GSE68857
g857 <- arrayDatasets[["GSE68857"]]
head(exprs(g857)) ## log=FALSE
g857.pd <- varpd(g857)
g857.pd$PRKCE <- with(g857.pd, grepl("PKC", characteristics_ch1.2))
g857.pd$Female <- with(g857.pd, grepl("female", characteristics_ch1.3))
g857.pd$Treatment <- with(g857.pd, factor(ifelse(grepl("FG-3149", characteristics_ch1.4), "IgG", "antiCTGF"),
                                          levels=c("IgG", "antiCTGF")))
g857.design <- with(g857.pd, model.matrix(~PRKCE*Female*Treatment))
colnames(g857.design) <- c("Baseline", "PRKCEtransgen", "Female", "antiCTGF", "PRKCE_Female", "PRKCE_antiCTGF", "Female_antiCTGF", "PRKCE_Female_antiCTGF")
g857.contrast <- makeContrasts("PRKCEtransgen",
                               levels=g857.design)
exportArray("GSE68857", log2=FALSE, design=g857.design, contrast=g857.contrast)

## GSE63759
g759 <- arrayDatasets[["GSE63759"]]
head(exprs(g759)) ## log=FALSE
g759.pd <- varpd(g759)
g759.design <- with(g759.pd, model.matrix(~source_name_ch1))
colnames(g759.design) <- c("Control","CoupTFII_OE")
g759.contrast <- makeContrasts("CoupTFII_OE", levels=g759.design)
exportArray("GSE63759", log2=FALSE, design=g759.design, contrast=g759.contrast)


## GSE67492
g492 <- arrayDatasets[["GSE67492"]]
head(exprs(g492)) ## log =FALSE
g492.pd <- varpd(g492)
(g492.design <- with(g492.pd, model.matrix(~characteristics_ch1)))
colnames(g492.design) <- c("control", "IDC", "BMPR2mut_PAH")
(g492.contrast <- makeContrasts("IDC", "BMPR2mut_PAH", levels=g492.design))
g492.features <- gsub(".*(ENST[0-9]*).*", "\\1", as.character(fData(g492)$gene_assignment))
g492.features[grepl("^-", g492.features)] <- NA
g492.isChar <- !is.na(g492.features) & !grepl("^ENST[0-9]*$", g492.features)
g492.features[g492.isChar] <- sapply(strsplit(g492.features[g492.isChar], "//"), function(x) trim(x[[1]]))
exportArray("GSE67492", log2=FALSE, design=g492.design, contrast=g492.contrast, featureID=g492.features)

## GSE63847
g847 <- arrayDatasets[["GSE63847"]]
head(exprs(GSE63847)) ## log=FALSE
g847.pd <- varpd(g847)
g847.design <- with(g847.pd, model.matrix(~characteristics_ch1*characteristics_ch1.4))
colnames(g847.design) <- c("Baseline", "C57BL.6J", "Trypomastigotes_150", "Trypomastigotes_200",
                           "C57BL.6J_T150", "C57BL.6J_T200")
g847.contrast <- makeContrasts("Trypomastigotes_150", "Trypomastigotes_200",
                               levels=g847.design)
exportArray("GSE63847", log2=FALSE, design=g847.design, contrast=g847.contrast, featureID=as.character(fData(g847)$GENE))

## GSE54681
g681 <- arrayDatasets[["GSE54681"]]
head(exprs(g681)) ## log=TRUE
g681.pd <- varpd(g681)
g681.pd$Group <- with(g681.pd, gsub("-[1-4]$", "", title))
g681.pd$Group[grepl("control", g681.pd$Group)] <- "control"
g681.pd$Group <- gsub("tamoxifen-day ([0-9]*)", "Tamoxifen_d\\1", g681.pd$Group)

g681.pd$Group <- factor(g681.pd$Group,
                        levels=c("control", "Tamoxifen_d2", "Tamoxifen_d3", "Tamoxifen_d10", "Tamoxifen_d28"))
(g681.design <- with(g681.pd, model.matrix(~Group)))
colnames(g681.design) <- levels(g681.pd$Group)
(g681.contrast <- makeContrasts(contrasts=levels(g681.pd$Group)[-1], levels=g681.design))
exportArray("GSE54681", log2=TRUE, design=g681.design, contrast=g681.contrast, featureID=fData(g681)$Entrez_Gene_ID)

## GSE82290
g290 <- arrayDatasets[["GSE82290"]]
head(exprs(g290)) ## log=TRUE
g290.pd <- varpd(g290)
g290.pd$isPatient <- with(g290.pd, !grepl("control", characteristics_ch1.1))
g290.pd$hasDCM <- with(g290.pd, grepl("dilated", characteristics_ch1.1) & isPatient) 
(g290.design <- with(g290.pd, model.matrix(~isPatient+hasDCM)))
colnames(g290.design) <- c("Baseline", "DCM", "conductionDefect")
(g290.contrast <- makeContrasts("DCM", levels=g290.design))
exportArray("GSE82290", log2=TRUE, design=g290.design, contrast=g290.contrast)

## GSE52601
g601 <- arrayDatasets[["GSE52601"]]
head(exprs(g601)) ## log=FALSE
(g601.pd <- varpd(g601))
g601.pd$Gender <- sapply(strsplit(as.character(g601$characteristics_ch1.3), ":"), function(x) trim(x[[2]]))
boxplot(exprs(g601)[grep("XIST", fData(g601)$Symbol),]~as.character(g601.pd$Gender)) ## NA is a female
g601.pd$Gender[g601.pd$Gender=="N/A"] <- "Female"
g601.pd$EF <- as.numeric(sapply(strsplit(as.character(g601$characteristics_ch1.6), ":"), function(x) trim(x[[2]])))
g601.pd$Disease <- sapply(strsplit(as.character(g601$characteristics_ch1.1), ":"), function(x) trim(x[[2]]))
g601.pd$Treatment <- sapply(strsplit(as.character(g601$characteristics_ch1.7), ":"), function(x) trim(x[[2]]))
g601.pd$Race <- sapply(strsplit(as.character(g601$characteristics_ch1.4), ":"), function(x) trim(x[[2]]))
g601.pd$Group <- factor(make.names(gsub("-rep[0-9]", "", g601$title)))

boxplot(g601.pd$EF~g601.pd$Group) ## EF is strongly correlated with the groups

g601.design <- with(g601.pd, model.matrix(~Group+Gender))
colnames(g601.design) <- c("Fetal", "DCM", "DCM.VAD", "ICM", "ICM.VAD", "nonfailing", "Male")
g601.contrast <- makeContrasts("DCM"="DCM-nonfailing",
                               "ICM"="ICM-nonfailing",
                               levels=g601.design)
exportArray("GSE52601", log2=FALSE, design=g601.design, contrast=g601.contrast)

## NGS

ngsFiles <- sapply(ngsIDs, function(x) dir(file.path(BASEDIR, x), pattern="*.txt$", full.name=TRUE))
ngsTbls <- lapply(ngsFiles, readMatrix)
names(ngsTbls) <- ngsIDs

exportNGS <- function(id, log2=TRUE, design, contrast, featureID=NULL) {
    stopifnot(id %in% names(ngsDs))
    eset <- ngsDs[[id]]
    exprs(eset) <- ngsTbls[[id]]
    if(log2 & any(exprs(eset)<=0)) {
        ## filter lowly expressed genes
        isLow <- apply(exprs(eset), 1, function(x) mean(x<=1)>=0.5 || any(x<=0))
        eset <- eset[!isLow,]
        if(!is.null(featureID))
            featureID <- featureID[!isLow]
    }
    exportStudy(eset, id, log2=log2, design=design, contrast=contrast, featureID=featureID)
}

## GSE71613
g613 <- ngsDs[["GSE71613"]]
g613.pd <- varpd(g613)
g613.pd$Group <- ofactor(gsub("diagnosis: ", "", g613.pd$ characteristics_ch1.1))
g613.design <- with(g613.pd, model.matrix(~Group))
colnames(g613.design) <- c("Healthy", "RestrictiveCM", "DCM")
(g613.contrast <- makeContrasts("RestrictiveCM", "DCM", levels=g613.design))
exportNGS("GSE71613", log2=TRUE, design=g613.design, contrast=g613.contrast)

## GSE71912
g912 <- ngsDs[["GSE71912"]]
g912.pd <- varpd(g912)
g912.pd$Group <- factor(with(g912.pd, ifelse(grepl("WT", characteristics_ch1.2), "WT", "Mib1_KO")),
                        levels=c("WT", "Mib1_KO"))
(g912.design <- with(g912.pd, model.matrix(~Group)))
colnames(g912.design) <- levels(g912.pd$Group)
(g912.contrast <- makeContrasts("Mib1_KO", levels=g912.pd$Group))
exportNGS("GSE71912", log2=FALSE, design=g912.design, contrast=g912.contrast)

## GSE65446
g446 <- ngsDs[["GSE65446"]]
summary(ngsTbls[["GSE65446"]])
(g446.pd <- varpd(g446))
(g446.design <- with(g446.pd, model.matrix(~source_name_ch1)))
colnames(g446.design) <- c("control", "hDCM")
(g446.contrast <- makeContrasts("hDCM", levels=g446.design))
exportNGS("GSE65446", log2=TRUE, design=g446.design, contrast=g446.contrast)

## GSE64391
g391 <- ngsDs[["GSE64391"]]
(g391.pd <- varpd(g391))
g391.pd$Group <- factor(with(g391.pd, ifelse(grepl("wildtype", source_name_ch1),
                                             "WT", "Bmi1_KO")),
                        levels=c("WT", "Bmi1_KO"))
(g391.design <- with(g391.pd, model.matrix(~Group)))
colnames(g391.design) <- c("control", "Bmi1_KO")
(g391.contrast <- makeContrasts("Bmi1_KO", levels=g391.design))
exportNGS("GSE64391", log2=TRUE, design=g391.design, contrast=g391.contrast)


##----------------------------------------##
## Perform DGE analysis with MPS signatures
##----------------------------------------##
getGSEid <- function(filename) sapply(strsplit(basename(filename), "-"), "[[", 1L)
DGEDATADIR <- file.path(BASEDIR, "data-for-DGE")
gseStudies <- getGSEid(dir(DGEDATADIR, pattern="*.gct"))
DGEOUTDIR <- file.path(BASEDIR, "DGE-outdirs")

buildLimmaComm <- function(gse) {
    exprsFile <- file.path(DGEDATADIR, sprintf("%s-expression.gct", gse))
    designFile <- file.path(DGEDATADIR, sprintf("%s-designMatrix.txt", gse))
    contrastFile <- file.path(DGEDATADIR, sprintf("%s-contrastMatrix.txt", gse))
    featureFile <- file.path(DGEDATADIR, sprintf("%s-featureAnnotation.txt", gse))
    outdir <- file.path(DGEOUTDIR, gse)
    logfile <- file.path(DGEOUTDIR, sprintf("%s-limma.log", gse))
    assertFile(c(exprsFile, designFile, contrastFile, featureFile))
    comm <- sprintf("ml load bi-R; /apps64/bi/R/R-3.3.0/bin/Rscript %s -infile %s -design %s -contrast %s -outdir %s -featureAnno %s -doPathwayAnalysis -mps -log %s",
                    "/apps64/bi/geneexpression/bin/maDge_limma.Rscript",
                    exprsFile,
                    designFile,
                    contrastFile,
                    outdir,
                    featureFile,
                    logfile)
    return(comm)
}

gseComms <- sapply(gseStudies, buildLimmaComm)

## run locally
## sapply(gseComms, system)
if(file.exists("limmaCommand.bash"))  file.remove("limmaCommand.bash")
submit <- sprintf("/apps64/bi/jobsubmitter/bin/submit_job -command \"%s\" -name %s -qsubFile limmaCommand.bash",
                  gseComms, gseStudies)
sapply(submit, system)
system(sprintf("ssh rbalhpc05 %s/limmaCommand.bash", getwd()))

## probmatic cases
#### featureAnno not used properly:GSE16909 GSE3586 GSE54681 GSE63847
#### design/contrast problem: GSE65446 GSE64391
## data.frame list into integer: GSE71613

