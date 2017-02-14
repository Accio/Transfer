## collect cardiomyopathy data for the PDD project

library(GEOquery)
library(limma)

library(ribiosUtils)
library(ribiosAnnotation)
library(ribiosIO)
library(ribiosPlot)
library(ribiosBioQC)

DATA_DIR <- "data"

f2i <- function(x) as.integer(as.character(x))
f2n <- function(x) as.numeric(as.character(x))
ofactor <- function(x) factor(x, levels=unique(x))
outfig <- function(x) file.path("figures", x)

cardioRData <- iofile("cardiomyopathy-meta.RData")
if(!loadFile(cardioRData)) {
  wittchen.raw <- getGEO("GSE4172")
  gaertner.raw <- getGEO("GSE29819")
  ameling.raw <- getGEO("GSE17800")
  ## DKFZ
  kittleson.raw <- getGEO("GSE1869")
  kuner.raw <- getGEO("GSE3586")
  barth.raw <- getGEO("GSE3585")
  
  harding.raw <- getGEO("GSE16909")
  witt.raw <- getGEO("GSE17749")
  glynjones.raw <- getGEO("GSE5606")

  li.raw <- getGEO("GSE51483")
  
  wittchen <- wittchen.raw[[1]]
  gaertner <- gaertner.raw[[1]]
  ameling <- ameling.raw[[1]]
  kittleson <- kittleson.raw[[1]]
  kuner <- kuner.raw[[1]]
  barth <- barth.raw[[1]]
  harding <- harding.raw[[1]]
  witt.1 <- witt.raw[[1]]
  witt.2 <- witt.raw[[2]]
  glynjones <- glynjones.raw[[1]]
  li <- li.raw[[1]]

  myNorm <- function(eset) {
    norm <- normalizeVSN(eset)
    exprs(eset) <- norm
    return(eset)
  }
  
  ## pheno data clean-ups
  ##----------------------------------------##
  ## Wittchen et al. 2007, DCMi data, Hans Lehrach group
  ## endomyocardial biopsy from the right ventricle, inflammation/PVB19 
  ##----------------------------------------##
  wittchen.pd <- data.frame(title=wittchen$title,
                            group=factor(ifelse(wittchen$source_name_ch1=="healthy control","healthy", "DCMi"),
                              levels=c("healthy", "DCMi")),
                            age=f2i(gsub("age: ", "", wittchen$characteristics_ch1)),
                            gender=gsub("gender: ", "", wittchen$characteristics_ch1.1),
                            EF=f2i(gsub("ejection fraction: ", "", wittchen$characteristics_ch1.2)),
                            LVEDD=f2i(gsub("left ventricular end diastolic diameter: ",
                              "", wittchen$characteristics_ch1.3)),
                            row.names=rownames(pData(wittchen)))
  wittchen.newPD <- new("AnnotatedDataFrame",
                        wittchen.pd,
                        varMetadata=data.frame(labelDescription=c("title","group","age","gender",
                          "eject fraction", "left ventricular end-diastolic diameter")))
  phenoData(wittchen) <- wittchen.newPD
  wittchen.fd <- annotateProbesets(featureNames(wittchen), chip="HG-U133_PLUS_2")
  fData(wittchen) <- wittchen.fd
  wittchen <- myNorm(wittchen)
  
  ##----------------------------------------##
  ## Gaertner et al., PG 2012
  ##----------------------------------------##
  gaertner.pd.backup <- pData(gaertner)
  gaertner.pd <- with(pData(gaertner),
                            data.frame(age=f2i(gsub("age: ", "", characteristics_ch1)),
                                       gender=gsub("gender: ", "", characteristics_ch1.2),
                                       ventricle=factor(gsub("ventricle: ", "", characteristics_ch1.3),levels=c("right", "left")),
                                       row.names=rownames(pData(gaertner))))
  gaertner.pd$gender[gaertner.pd$gender=="n.a."] <- NA
  gaertner.pd$gender[gaertner.pd$gender=="w"] <- "f"
  gaertner.pd$gender <- factor(as.character(gaertner.pd$gender), levels=c("f", "m"))
  levels(gaertner.pd$gender) <- c("female", "male")
  gaertner.pd$group <- factor(sapply(strsplit(as.character(gaertner$source_name_ch1), " "), "[[", 4),
                              levels=c("NF", "DCM", "ARVC"))
  gaertner.pd$patient <- paste(gaertner.pd$group,
                               sapply(strsplit(as.character(gaertner$source_name_ch1), " "), "[[", 6),
                               sep=".")
  gaertner.newPD <- new("AnnotatedDataFrame",
                        gaertner.pd,
                        varMetadata=data.frame(labelDescription=c("Age", "Gender", "ventricle",
                                                 "group: ARVC=Arrhythmogenic right ventricular cardiomyopathy; DCM=idiopathic dilated cardiomyopathy; NF=not failing",
                                                 "patient: unique patient identifiers")))
  phenoData(gaertner) <- gaertner.newPD
  gaertner.fd <- annotateProbesets(featureNames(gaertner), "HG-U133_PLUS_2")
  fData(gaertner) <- gaertner.fd
  gaertner <- myNorm(gaertner)
  ## exprs(gaertner)[grep("XIST", fData(gaertner)$GeneSymbol),34:38]
  gaertner$gender[35:36] <- c("male", "male")
  gaertner$age[35:36] <- mean(gaertner$age, na.rm=TRUE)
  
  ##----------------------------------------##
  ## Ameling et al., EHJ 2012
  ## retrospective analysis of IA/IgG (Immunoadsorption with subsequent immunoglobulin G substitution) responders/non-responders
  ##----------------------------------------##
  ameling.group <- rep(NA, nrow(pData(ameling)))
  ameling.group[grepl("Co", ameling$source_name_ch1)] <- "control"
  ameling.pd <- with(pData(ameling),
                     data.frame(LVEF=f2i(gsub("lvef: ", "", characteristics_ch1.1)),
                                LVIDd=f2i(gsub("lvidd: ", "", characteristics_ch1.2)),
                                gender=gsub("gender: ","", characteristics_ch1.3),
                                age=f2i(gsub("age: ", "", characteristics_ch1.4)),
                                virus=gsub("virus: ", "",  characteristics_ch1.5),
                                inflammation=f2i(gsub("inflammation \\(cd68\\+ \\+ cd3\\+\\): ", "", characteristics_ch1.6)),
                                bmi=f2i(gsub("bmi: ", "", characteristics_ch1.7))))
  ameling.fd <- annotateProbesets(featureNames(ameling), chip="HG-U133_PLUS_2")
  ameling <- myNorm(ameling)
  fData(ameling) <- ameling.fd
  ## NOT FINISHED: responder/non-responder, and disease duration is missing

  ##----------------------------------------##
  ## Kittleson et al., PG 2005
  ## ischemic and nonischemic cardiomyopathy
  ##----------------------------------------##
  kittleson.pdbackup <- pData(kittleson)
  kittleson.group <- rep(NA, nrow(kittleson.pdbackup))
  kittleson.group[kittleson$description=="Explanted heart, nonischemic cardiomyopathy"] <- "NICM"
  kittleson.group[kittleson$description=="Expanted heart, nonischemic cardiomyopathy"] <- "NICM"
  kittleson.group[kittleson$description=="Expanted heart, ischemic cardiomyopathy"] <- "ICM"
  kittleson.group[kittleson$description=="Pre-LVAD, ischemic cardiomyopathy"] <- "PreLVAD.ICM"
  kittleson.group[kittleson$description=="Pre-LVAD, nonischemic cardiomyopathy"] <- "PreLVAD.NICM"
  kittleson.group[kittleson$description=="Unused donor heart"] <- "Nonfailing"
  kittleson.group <- factor(kittleson.group,
                            levels=c("Nonfailing", "ICM", "NICM", "PreLVAD.ICM", "PreLVAD.NICM"))
  
  kittleson.pd <- data.frame(title=kittleson$title, group=kittleson.group,
                             LVAD=grepl("LVAD", kittleson.group),
                             NICM=grepl("NICM", kittleson.group),
                             row.names=sampleNames(kittleson))
  kittleson.newPD <- new("AnnotatedDataFrame",
                         kittleson.pd,
                         varMetadata=data.frame(labelDescription=c("Title",
                                                  "Group: ICM=ischemic cardiomyopathy;NCIM:nonischemic cardiomyopathy;VLAD:left ventricular assist device",
                                                  "pre-LVAD", "nonischemic cardiomyopathy")))
  phenoData(kittleson) <- kittleson.newPD
  kittleson.fd <- annotateProbesets(featureNames(kittleson), "HG-U133_PLUS_2")
  fData(kittleson) <- kittleson.fd
  
  ##----------------------------------------##
  ## Barth et al., JACC, 2006 (kuner)
  ## non-failing versus DCM
  ##----------------------------------------##
  kuner.pdbackup <- pData(kuner)
  kuner.group <- factor(ifelse(kuner$characteristics_ch1=="Non failing", "Nonfailing", "DCM"),
                        levels=c("Nonfailing", "DCM"))
  kuner.pd <- data.frame(group=kuner.group,
                         row.names=rownames(pData(kuner)))
  kuner.newPD <- new("AnnotatedDataFrame",
                     kuner.pd)
  phenoData(kuner) <- kuner.newPD
  annotateBarth <- function(fd) {
    refseqs <- sapply(strsplit(as.character(fd$GB_LIST), " "), function(x) {
      if(length(x)>=1)
        return(x[1])
      return(NA)
    })
    anno0 <- annotatemRNAs(refseqs)
    gs <- annotateGeneSymbols(as.character(fd$GENE_SYMBOL))
    isFromGs <- is.na(anno0$GeneID) & !is.na(gs$GeneID)
    geneids <- anno0$GeneID; geneids[isFromGs] <- gs$GeneID[isFromGs]
    genesymbols <- anno0$GeneSymbol; genesymbols[isFromGs] <- gs$GeneSymbol[isFromGs]
    genename <- anno0$GeneName; genename[isFromGs] <- gs$GeneName[isFromGs]
    res <- data.frame(mRNA=fd$REFSEQ_ID,
                      RawGeneSymbol=fd$GENE_SYMBOL,
                      GeneID=geneids,
                      GeneSymbol=genesymbols,
                      GeneName=genename)
    rownames(res) <- rownames(fd)
    return(res)
  }
  kuner.fd <- annotateBarth(fData(kuner))
  fData(kuner) <- kuner.fd
  
  ##----------------------------------------##
  ## Barth et al., JACC, 2006 (Barth)
  ## non-failing versus DCM
  ##----------------------------------------##
  barth.pdbackup <- pData(barth)
  barth.group <- factor(ifelse(barth$characteristics_ch1=="Non failing", "Nonfailing", "DCM"),
                        levels=c("Nonfailing", "DCM"))
  barth.pd <- data.frame(group=barth.group,
                         row.names=rownames(pData(barth)))
  barth.newPD <- new("AnnotatedDataFrame",
                     barth.pd)
  phenoData(barth) <- barth.newPD
  barth.fd <- annotateProbesets(featureNames(barth), "HG-U133_PLUS_2")
  fData(barth) <- barth.fd
  
  ##----------------------------------------##
  ## Harding et al.
  ## EP4 KO mice model 
  ##----------------------------------------##
  harding.pdbackup <- pData(harding)
  harding.age <- f2n(gsub("age \\(weeks\\): ", "", harding$characteristics_ch1.1))
  (harding.ef <- f2n(gsub(".*: ", "", harding$characteristics_ch1.2)))
  harding.strain <- factor(ifelse(grepl("WT", harding$title), "WT", "EP4_KO"),
                           levels=c("WT", "EP4_KO"))
  harding.group <- as.character(harding.strain)
  harding.group[grepl("low EF", harding$title)] <- "EP4_KO_LowEF"
  harding.group[grepl("high EF", harding$title)] <- "EP4_KO_HighEF"
  harding.pd <- data.frame(Age=harding.age,
                           EF=harding.ef,
                           strain=harding.strain,
                           group=harding.group,
                           GEOdescription=harding$description)
  harding.newPD <- new("AnnotatedDataFrame",
                       harding.pd,
                       varMetadata=data.frame(labelDescription=c("Age","Ejection Fraction",
                                                "Strain", "Group", "GEO description (not knowing whether its relevant or not")))
  phenoData(harding) <- harding.newPD

  ##----------------------------------------##
  ## Vasti et al. [witt]
  ##----------------------------------------##
  witt.1.pdbackup <- pData(witt.1)
  witt.1.group1 <- ifelse(grepl("Doxorubicin", witt.1$source_name_ch1), "Doxorubicin", "control")
  witt.1.group2 <- ifelse(grepl("Doxorubicin", witt.1$source_name_ch2), "Doxorubicin", "control")
  witt.1.newPD <- new("AnnotatedDataFrame",
                      data.frame(group.ch1=witt.1.group1, group.chr2=witt.1.group2, row.names=sampleNames(witt.1)))
  phenoData(witt.1) <- witt.1.newPD

  witt.2.pdbackup <- pData(witt.2)
  witt.2.group1 <-  c("ErbB4_KO", "WT", "ErbB4_KO", "WT", "ErbB4_KO_Doxorubicin","WT", "ErbB4_KO_Doxorubicin", "WT")
  witt.2.group2 <-  c("WT", "ErbB4_KO", "WT", "ErbB4_KO", "WT", "ErbB4_KO_Doxorubicin","WT", "ErbB4_KO_Doxorubicin")
  witt.2.newPD <- new("AnnotatedDataFrame",
                      data.frame(group.ch1=witt.2.group1, group.chr2=witt.2.group2, row.names=sampleNames(witt.2)))
  phenoData(witt.2) <- witt.2.newPD

  ##----------------------------------------##
  ## Glyn-Jones et al., PG 2006
  ## Rat diabetic cardiomyopathy model
  ##----------------------------------------##
  glynjones.pdbackup <- pData(glynjones)
  glynjones.newPD <- new("AnnotatedDataFrame",
                         data.frame(group=factor(ifelse(grepl("Diabetic", glynjones$source_name_ch1), "Diabetic", "control"),
                                      levels=c("control", "Diabetic")),
                                    row.names=sampleNames(glynjones)))
  phenoData(glynjones) <- glynjones.newPD
  glynjones.fd <- annotateProbesets(featureNames(glynjones), "RAT230_2", orthologue=TRUE)
  fData(glynjones) <- glynjones.fd
  glynjones <- myNorm(glynjones)

  ##----------------------------------------##
  ## Li et al., PG 2014
  ## mouse stem-cell model
  ##----------------------------------------##
  li.pdbackup <- pData(li)
  li.type <- sapply(strsplit(as.character(li$characteristics_ch1), ":"), "[[", 1L)
  li.sample <- li$characteristics_ch1
  li.sample[li.sample=="tissue: Right ventricle tissue"] <- "tissue: Right ventricle tissues"
  li.sample[li.sample=="tissue: Left ventricle tissue"] <- "tissue: Left ventricle tissues"
  li.sample <- ofactor(li.sample)
  levels(li.sample) <- c("cellLineR1", "wholeEmbryoTissue", "wholeHeartTube", "leftVentricle", "rightVentricle")
  li.age <- ofactor(gsub("age: ","", li$characteristics_ch1.1))
  levels(li.age) <- c("NA", "E7.5", "E8.5", "E9.5", "E12.5", "E14.5", "E18.5", "D3AfterBirth", "adult")
  li.group <- ofactor(li$description)
  li.pdNew <- new("AnnotatedDataFrame",
                  data.frame(group=li.group,
                             type=li.type,
                             sample=li.sample,
                             age=li.age,
                             row.names=sampleNames(li)))
  phenoData(li) <- li.pdNew
  li.fd <- annotateProbesets(featureNames(li), "MOUSE430_2", orthologue=TRUE)
  fData(li) <- li.fd

  ##----------------------------------------##
  ## statistical models
  ##----------------------------------------##
  wittchen.model <- model.matrix(~group+age+gender+EF+LVEDD, data=pData(wittchen))
  colnames(wittchen.model) <- c("baseline", "DCMi", "age", "gender", "EF", "LVEDD")
  wittchen.contrast <- makeContrasts(DCMi, EF, LVEDD, levels=colnames(wittchen.model))

  #### very difficult to take block factor into consideration
  ## gaertner.nind <- ncol(gaertner)/2
  ## gaertner.ind <- gl(gaertner.nind,2)
  ## gaertner.model <- model.matrix(~gaertner.ind+ventricle, data=pData(gaertner))
  ## lmFit(gaertner, gaertner.model)
  
  gaertner.model <- model.matrix(~group+gender+ventricle+group:ventricle+age, data=pData(gaertner))
  colnames(gaertner.model) <- c("baseline", "DCM", "ARVC", "male", "leftVentricle", "age", "DCM_leftVentricle", "ARVC_leftVentricle")
  gaertner.contrast <- makeContrasts(DCM_right="DCM", ARVC_right="ARVC",
                                     NF_left="leftVentricle",
                                     DCM_left=DCM+DCM_leftVentricle, ARVC_left=ARVC+ARVC_leftVentricle,
                                     levels=colnames(gaertner.model))
  
  kittleson.model <- model.matrix(~LVAD*(!NICM), data=pData(kittleson))
  colnames(kittleson.model) <- c("baseline", "LVAD", "ICM", "LVAD_and_ICM")
  kittleson.contrast <- makeContrasts(ICM, LVAD_and_ICM, LVAD, levels=colnames(kittleson.model))
  
  kuner.model <- model.matrix(~group, data=pData(kuner))
  colnames(kuner.model) <- c("baseline", "DCM")
  kuner.contrast <- makeContrasts(DCM, levels=colnames(kuner.model))
  
  barth.model <- model.matrix(~group, data=pData(barth))
  colnames(barth.model) <- c("baseline", "DCM")
  barth.contrast <- makeContrasts(DCM, levels=colnames(barth.model))
  
  harding.model <- model.matrix(~strain+EF+Age, data=pData(harding))
  colnames(harding.model) <- c("baseline", "EP4_KO", "EF", "Age")
  harding.contrast <- makeContrasts(EP4_KO, EF, levels=colnames(harding.model))
  
  glynjones.model <- model.matrix(~group, data=pData(glynjones))
  colnames(glynjones.model) <- c("baseline", "Diabetic")
  glynjones.contrast <- makeContrasts(Diabetic, levels=colnames(glynjones.model))


  ## li's model is complicated
  li.isTissue <- as.integer(grepl("tissue", li$group))
  li.isLeft <- as.integer(li$sample=="leftVentricle")
  li.isRight <- as.integer(li$sample=="rightVentricle")
  ## li.wholeembyro <- as.integer(grepl("whole embryo", li$group))
  ## li.wholeheart <- as.integer(grepl("whole heart", li$group, ignore.case=TRUE))
  li.E7.5 <- as.integer(li$age=="E7.5") ## whole embryo
  li.E8.5 <- as.integer(li$age=="E8.5") ## whole heart
  li.E9.5 <- as.integer(li$age=="E9.5")
  li.E12.5 <- as.integer(li$age=="E12.5")
  li.E14.5 <- as.integer(li$age=="E14.5")
  li.E18.5 <- as.integer(li$age=="E18.5")
  li.d3 <- as.integer(li$age=="D3AfterBirth")
  li.adult <- as.integer(li$age=="adult")
  li.model <- model.matrix(~li.isTissue+li.E8.5+
                           li.E9.5+li.E12.5+li.E14.5+li.E18.5+li.d3+li.adult+
                           (li.E9.5+li.E12.5+li.E14.5+li.E18.5+li.d3+li.adult):li.isRight)
  colnames(li.model) <- c("baseline", "tissue", "E8.5",
                          "E9.5", "E12.5", "E14.5", "E18.5", "d3", "adult", ## everything compared to E7.5
                          "E9.5right", "E12.5right", "E14.5right", "E18.5right", "d3right", "adultright")
  poly6 <- contr.poly(6)
  li.contrast.raw <- makeContrasts(tissue=tissue,
                               E8.5=E8.5,
                               E9.5.vs.E8.5=E9.5-E8.5,
                               E12.5.vs.E9.5=E12.5-E9.5,
                               E14.5.vs.E12.5=E14.5-E12.5,
                               E18.5.vs.E14.5=E18.5-E14.5,
                               d3.vs.E18.5=d3-E18.5,
                               adult.vs.d3=adult-d3,
                               E9.5right.vs.E8.5=E9.5right-E8.5,
                               E12.5right.vs.E9.5right=E12.5right-E9.5right,
                               E14.5right.vs.E12.5right=E14.5right-E12.5right,
                               E18.5right.vs.E14.5right=E18.5right-E14.5right,
                               d3right.vs.E18.5right=d3right-E18.5right,
                               adultright.vs.d3right=adultright-d3right, levels=colnames(li.model))
  li.contrast.poly <- cbind(leftL=c(rep(0,3),poly6[, 1L], rep(0, 6)),
                            leftQ=c(rep(0,3),poly6[, 2L], rep(0, 6)),
                            rightL=c(rep(0,3), poly6[, 1L], poly6[, 1L]),
                            rightQ=c(rep(0,3), poly6[, 2L], poly6[, 2L]))
  li.contrast <- cbind(li.contrast.raw, li.contrast.poly)
 ##li.ageMat.raw <- model.matrix(~li.isRight:age, pData(li))
 ##
 ##li.model <- cbind(baseline=1, li.isCellLine, li.isLeft, li.isRight, li.ageMat.raw[,-1])
 ##li.model <- li.model.raw[,apply(li.model.raw, 2, function(x) any(x!=0))]
 ##colnames(li.model) <- c("baseline", "E7.5", "E8.5", "E9.5", "E12.5", "E14.5", "E18.5", "D3AfterBirth", "Adult",
 ##                        "wholeEmbryoTissue", "WholeHeartTube", "leftVentricle", "rightVentricle",
 ##                        "wholeEmbryoTissue_E7.5", "
                          
  ##----------------------------------------##
  ## Make fits
  ##----------------------------------------##

  myfit <- function(eset, design, contrast) {
    mat <- exprs(eset)
    weights <- arrayWeights(eset, design)
    fit <- lmFit(eset, design, weights=weights)
    fit2 <- contrasts.fit(fit, contrasts=contrast)
    efit <- eBayes(fit2)
    efit$genes <- fData(eset)
    topTables <- lapply(1:ncol(contrast), function(i) {
      topTable(efit, coef=i, number=nrow(eset), sort.by="p")
    })
    colnames(topTables) <- colnames(contrast)
    res <- list(weights=weights,
                efit=efit,
                topTables=topTables)
  }
  
  
  wittchen.fit <- myfit(wittchen, wittchen.model, wittchen.contrast)
  gaertner.fit <- myfit(gaertner, gaertner.model, gaertner.contrast)
  ## ameling was left out due to missing parameters
  kittleson.fit <- myfit(kittleson, kittleson.model, kittleson.contrast)
  kuner.fit <- myfit(kuner, kuner.model, kuner.contrast)
  barth.fit <- myfit(barth, barth.model, barth.contrast)
  harding.fit <- myfit(harding, harding.model, harding.contrast)
  ## witt.1 and witt.2 left out due to two-channel design
  glynjones.fit <- myfit(glynjones, glynjones.model, glynjones.contrast)
  li.fit <- myfit(li, li.model, li.contrast)
  
  ##----------------------------------------##
  ## save the datasets
  ##----------------------------------------##
  save(wittchen, gaertner, ameling,
       kittleson, kuner, barth,
       harding, witt.1, witt.2, glynjones,
       li,
       
       wittchen.model, wittchen.contrast,
       gaertner.model, gaertner.contrast,
       kittleson.model, kittleson.contrast,
       kuner.model, kuner.contrast,
       barth.model, barth.contrast,
       harding.model, harding.contrast,
       glynjones.model, glynjones.contrast,
       li.model, li.contrast,

       wittchen.fit, gaertner.fit, kittleson.fit,
       kuner.fit, barth.fit, harding.fit, glynjones.fit, li.fit,
       file=cardioRData)
}


## gene set enrichment analysis
tt.kittleson.ICM <- kittleson.fit$topTables[[1L]]
tt.kuner.DCM <- kuner.fit$topTables[[1L]]
tt.barth.DCM <- barth.fit$topTables[[1L]]

tt.wittchen.DCMi <- wittchen.fit$topTables[[1]]
## tt.wittchen.EF <- wittchen.fit$topTables[[2]]
## tt.wittchen.LVEDD <- wittchen.fit$topTables[[3]]
tt.gaertner.DCMright <- gaertner.fit$topTables[[1]]
tt.gaertner.DCMleft <- gaertner.fit$topTables[[4]]
write.table(tt.gaertner.DCMleft, outfig("topTable_gaertner_DCMleft.txt"), quote=FALSE, sep="\t")
write.table(tt.gaertner.DCMright, outfig("topTable_gaertner_DCMright.txt"), quote=FALSE, sep="\t")

tt.harding.EP4_KO <- harding.fit$topTables[[1]]

tt.li.leftL <- li.fit$topTables[[match("leftL", colnames(li.contrast))]]
tt.li.rightL <- li.fit$topTables[[match("rightL", colnames(li.contrast))]]

tt.glynjones.diabetic <- glynjones.fit$topTables[[1]]
matchPlot <- function(tbl1, tbl2,...) {
  probes <- unique(tbl1[, "ProbeID"], tbl2[, "ProbeID"])
  mt1 <- matchColumn(probes, tbl1, "ProbeID")
  mt2 <- matchColumn(probes, tbl2, "ProbeID")
  plot(mt1$logFC, mt2$logFC, ...)
  legend("topleft",
       sprintf("cor=%.2f",
               cor(mt1$logFC, mt2$logFC)), bty="n")
  abline(h=0, v=0)
}
matchPlot(tt.gaertner.DCMleft,
          tt.gaertner.DCMright,
          xlab="logFC [left ventricle]", ylab="logFC [right ventricle]", main="Gaertner et al.")
ipdf(outfig("Gaertner.left.right.xyplot.pdf"), width=6L, height=6L)
matchPlot(tt.li.leftL,
          tt.li.rightL,
          xlab="logFC [left ventricle]", ylab="logFC [right ventricle]", main="Gaertner et al.")
ipdf(outfig("Li.left.right.xyplot.pdf"), width=6L, height=6L)

myScatter <- function(tbl, ylim,  xlim, xdesc, ...) {
  isSig <- tbl$adj.P.Val<0.05
  ty <- abs(log10(tbl$P.Value))
  tx <- tbl$logFC
  if(missing(ylim))
    ylim <- c(0, pmax(2, max(ty)))
  if(missing(xlim))
    xlim <- c(-abs(max(tx)), abs(max(tx)))
  if(missing(xdesc)) {
    xlab <- "LogFC"
  } else {
    xlab <- sprintf("LogFC (%s)", xdesc)
  }
  smoothScatter(ty~tx, xlim=xlim, ylim=ylim,
                xlab=xlab, ylab="-log10(P)", ...)
  points(tx[isSig], ty[isSig], col="red", pch=1, cex=0.75)
}

op <- par(mfrow=c(3,2), mgp=c(2.5,1,0), mar=c(4,4,2,1)+0.5, cex.axis=1.25, cex.lab=1.5, cex.main=1.5)
myScatter(tt.wittchen.DCMi, main="Wittchen et al.", xdesc="DCMi vs ctrl")
myScatter(tt.kittleson.ICM, main="Kittleson et al.", xdesc="ICM vs ctrl")
myScatter(tt.kuner.DCM, main="Barth et al. [Kuner]", xdesc="DCM vs ctrl")
myScatter(tt.barth.DCM, main="Barth et al. [Barth]", xdesc="DCM vs ctrl")
myScatter(tt.gaertner.DCMleft, main="Gaertner et al.[left]", xdesc="DCM vs ctrl, left ventricle")
myScatter(tt.gaertner.DCMright, main="Gaertner et al.[right]", xdesc="DCM vs ctrl, right ventricle")

ipdf(outfig("DCM_humanStudies.pdf"), width=8, height=10)

par(mfrow=c(2,2), mgp=c(2.5,1,0), mar=c(4,4,2,1)+0.5, cex.axis=1, cex.lab=1.25, cex.main=1.25)
myScatter(tt.harding.EP4_KO, main="Harding et al.", xdesc="EP4-KO vs ctrl")
myScatter(tt.li.leftL, main="Li et al. [left]",xdesc="Linear changes, left ventricle")
myScatter(tt.li.rightL, main="Li et al. [right ventricle]", xdesc="Linear changes, right ventricle")
myScatter(tt.glynjones.diabetic, main="Glyn-Jones et al.", xdesc="Diabetic vs ctrl")
ipdf(outfig("DCM_animalModels.pdf"), width=7, height=7)

## human gene changes testedin the six studies
common.hg <- mintersect(fData(kittleson)$GeneID,
                        fData(kuner)$GeneID,
                        fData(barth)$GeneID,
                        fData(wittchen)$GeneID,
                        fData(gaertner)$GeneID)
length(common.hg) ## 5623 genes
lfc.hg <- cbind(Kittleson=matchColumn(common.hg, tt.kittleson.ICM, "GeneID")$logFC,
                Kuner=matchColumn(common.hg, tt.kuner.DCM, "GeneID")$logFC,
                Barth=matchColumn(common.hg, tt.barth.DCM, "GeneID")$logFC,
                Wittchen=matchColumn(common.hg, tt.wittchen.DCMi, "GeneID")$logFC,
                Gaertner.left=matchColumn(common.hg, tt.gaertner.DCMleft, "GeneID")$logFC,
                Gaertner.right=matchColumn(common.hg, tt.gaertner.DCMright, "GeneID")$logFC)
rownames(lfc.hg) <- common.hg

## PCA of changes
write.table(lfc.hg, outfig("lfc_table.txt"))
par(mfrow=c(1,1), mar=c(4,4,2,1), mgp=c(2,1,0))
plotPCA(prcomp(t(lfc.hg), scale=TRUE, center=TRUE), text=TRUE, points=FALSE,
        xlim=c(-0.3, 0.9), main="PCA of gene expression changes")
ipdf(outfig("PCA-gene-logFC.pdf"), width=5L, height=5L)
plotPCA(prcomp(t(lfc.hg[,-match("Wittchen", colnames(lfc.hg))]), scale=TRUE, center=TRUE),
        text=TRUE, points=FALSE, xlim=c(-0.5, 0.75))
ipdf(outfig("PCA-gene-logFC-excludingWittchen.pdf"), width=6L, height=6L)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
mypanel <- function(x,y,...) {
  abline(h=0, v=0, col="gray")
  panel.smooth(x,y,...)
}
pairs(lfc.hg[,c("Kittleson", "Kuner", "Barth", "Gaertner.left", "Gaertner.right", "Wittchen")],
      lower.panel = mypanel, upper.panel = panel.cor)
ipdf(outfig("gene-logFC-pairs.pdf"), width=8L, height=8L)
par(mar=c(3,3,2,1), mgp=c(2,1,0))
plot(hclust(dist(t(lfc.hg)), method="ward.D"),xlab="Euclidian distances between logFCs")
ipdf(outfig("gene-logFC-hclust.pdf"), width=4L, height=4L)


## overlap of common.hg with signature genes
sigGenes <- readLines("/homebasel/biocomp/zhangj83/projects/2012-11-RepPath/version2-2014/version2-genesymbols")
sigIDs <- annotateGeneSymbols(sigGenes)$GeneID
compTwoVecs(sigIDs, common.hg)
sig.lfc.hg <- lfc.hg[na.omit(intersect(sigIDs, common.hg)),]
plotPCA(prcomp(t(sig.lfc.hg)), text=TRUE, points=FALSE,
        xlim=c(-0.4, 0.9), main="PCA of signature gene expression changes")
ipdf(outfig("PCA-signatureGenes-logFC.pdf"), width=5L, height=5L)


## Gaertner's data
ronet.file <- "/data/bi/httpd_8080/htdoc/apps/gsea/genesets/path.ronet.symbol.roche.gmt"
upstream.file <- "/data/bi/httpd_8080/htdoc/apps/gsea/genesets/metabase.down.expression.symbol.roche.gmt"

ronet <- read_gmt_list(ronet.file)
upstream <- read_gmt_list(upstream.file)

matchGenes <- function(gmtlist, genesymbols, min=5L, max=2000L) {
  res <- lapply(gmtlist, function(x) unique(na.omit(match(x$genes, genesymbols))))
  res.len <- sapply(res, length)
  res <- res[res.len>=min & res.len<=max]
  return(res)
}
ronet.ind <- matchGenes(ronet, fData(gaertner)$GeneSymbol)
upstream.ind <- matchGenes(upstream, fData(gaertner)$GeneSymbol)

camera.ronet.NFleft <- camera(gaertner, index=ronet.ind,
                                design=gaertner.model,
                                contrast=gaertner.contrast[,"NF_left"])
camera.ronet.DCMright <- camera(gaertner, index=ronet.ind,
                                design=gaertner.model,
                                contrast=gaertner.contrast[,"DCM_right"])
camera.ronet.DCMleft <- camera(gaertner, index=ronet.ind,
                               design=gaertner.model,
                               contrast=gaertner.contrast[,"DCM_left"])
matchPathPlot <- function(tbl1, tbl2,...) {
  probes <- unique(rownames(tbl1), rownames(tbl2))
  mt1 <- matchColumn(probes, tbl1, 0L)
  mt2 <- matchColumn(probes, tbl2, 0L)
  mt1$score <- -log10(mt1$PValue) * ifelse(mt1$Direction=="Down", -1, 1)
  mt2$score <- -log10(mt2$PValue) * ifelse(mt2$Direction=="Down", -1, 1)
  plot(mt1$score, mt2$score, ...)
  legend("topleft",
       sprintf("cor=%.2f",
               cor(mt1$score, mt2$score, use="complete.obs")), bty="n")
  abline(h=0, v=0)
  abline(0, 1)
}

camera.upstream.DCMright <- camera(gaertner, index=upstream.ind,
                                design=gaertner.model,
                                contrast=gaertner.contrast[,"DCM_right"])
camera.upstream.DCMleft <- camera(gaertner, index=upstream.ind,
                                  design=gaertner.model,
                                  contrast=gaertner.contrast[,"DCM_left"])

par(mar=c(3,3,1,1)+1, mgp=c(2,1,0))
matchPathPlot(camera.ronet.DCMright, camera.ronet.DCMleft,
              xlab="Enrichment score [DCM right ventricle]",
              ylab="Enrichment score [DCM left ventricle]",
              main="RONET pathways")
ipdf(outfig("gaertner-ronet-scatter.pdf"), width=5L, height=5L)
matchPathPlot(camera.upstream.DCMright, camera.upstream.DCMleft,
              xlab="Enrichment score [DCM right ventricle]",
              ylab="Enrichment score [DCM left ventricle]",
              main="Upstream regulators")
ipdf(outfig("gaertner-upstream-scatter.pdf"), width=5L, height=5L)

## camera.summary
gaertner.camera <- rbind(cbind(type="RONET", sample="right", geneset=rownames(camera.ronet.DCMright), camera.ronet.DCMright),
                         cbind(type="RONET", sample="left", geneset=rownames(camera.ronet.DCMleft),camera.ronet.DCMleft),
                         cbind(type="Upstream", sample="right", geneset=rownames(camera.upstream.DCMright), camera.upstream.DCMright),
                         cbind(type="Upstream", sample="left", geneset=rownames(camera.upstream.DCMleft),camera.upstream.DCMleft))
gaertner.camera$ES <- abs(log10(gaertner.camera$PValue)) * ifelse(gaertner.camera$Direction=="Up", 1L, -1L)
rownames(gaertner.camera) <- NULL

gaertner.tbl <- reshape(gaertner.camera, idvar=c("type", "geneset", "NGenes"), timevar="sample", direction="wide")
xyplot(ES.left~ES.right | type, gaertner.tbl,
       scales=list(tck=c(1,0), alternating=1), layout=c(1,2),
       panel=function(x,y,...) {
         panel.xyplot(x,y,...)
         panel.abline(v=0, h=0)
         panel.abline(0,1)
       }, xlab="Enrichment score [right ventricle]", ylab="Enrichment score [left ventricle]")
ipdf(outfig("gaertner-scatters.pdf"), width=4L, height=7L)

write.table(gaertner.tbl, outfig("gaertner-ES-table.txt"))
gaertner.sigtbl <- subset(gaertner.tbl, abs(ES.left)>=3 | abs(ES.right)>=3)


## PLAGE
ronet.list <- lapply(ronet.ind, function(x) featureNames(gaertner)[x])
upstream.list <- lapply(upstream.ind, function(x) featureNames(gaertner)[x])
library(GSVA)
plage.ronet.DCM <- gsva(expr=gaertner, gset.idx.list=ronet.list,  method="plage")
plage.upstream.DCM <- gsva(expr=gaertner, gset.idx.list=upstream.list,  method="plage")

fitPLAGE <- function(plage, model, contrasts) {
  eset <- plage$es.obs
  fit <- lmFit(eset, model)
  fit <- contrasts.fit(fit, contrasts=contrasts)
  fit <- eBayes(fit)

  tbls <- do.call(rbind,lapply(1:ncol(contrasts), function(x) {
    res <- topTable(fit, coef=x, number=nrow(eset))
    res$pathway <- rownames(res)
    return(res)
  }))
  tbls$contrast <- rep(colnames(contrasts), each=nrow(eset))
  return(tbls)
}
plage.ronet.tables <- fitPLAGE(plage.ronet.DCM, gaertner.model, gaertner.contrast)
plage.upstream.tables <- fitPLAGE(plage.upstream.DCM, gaertner.model, gaertner.contrast)
plage.results <- rbind(cbind(type="RONET", plage.ronet.tables),
                       cbind(type="Upstream", plage.upstream.tables))
write.table(plage.results, outfig("gaertner-PLAGE-table.txt"))

hedgeInd <- ronet.ind[["Signaling_events_mediated_by_the_Hedgehog_family_NCI_NATURE"]]
nfkb <- ronet.ind[["Canonical_NF_kappaB_pathway_NCI_NATURE"]]
p53 <- ronet.ind[["p53_Dependent_G1_DNA_Damage_Response_REACTOME"]]
wnt <- ronet.ind[["deactivation_of_the_beta_catenin_transactivating_complex_REACTOME"]]
fn1 <- upstream.ind[["FN1_Influence on expression_Activation_Downstream"]]
acyl <- ronet.ind[["Acyl_chain_remodelling_of_PC_REACTOME"]]
g1s <- ronet.ind[["Cyclin_E_associated_events_during_G1_S_transition_REACTOME"]]

gaertner.left <- gaertner[,gaertner$ventricle=="left" & gaertner$group!="ARVC"]
gaertner.right <- gaertner[,gaertner$ventricle=="right" & gaertner$group!="ARVC"]
myHeatmap <- function(eset,...) {
  biosHeatmap(exprs(eset), zlim=c(-2,2),
              cexCol=1.5, cexRow=1.25,
              labRow=fData(eset)$GeneSymbol,
              scale="row", Colv=TRUE, dendrogram="row",
              color.key.title="Rel.Exp.", lhei=c(1,6),
              ColSideColor=c("yellow", "red", "black")[as.integer(eset$group)],
              ...)
}
pdf(outfig("pathway-heatmaps.pdf"), width=4L, height=6L)
myHeatmap(gaertner.left[hedgeInd,], main="Hedgehog/left")
myHeatmap(gaertner.right[hedgeInd,], main="Hedgehog/right")
myHeatmap(gaertner.left[nfkb,], main="NFkB/left")
myHeatmap(gaertner.right[nfkb,], main="NFkB/right")
myHeatmap(gaertner.left[p53,], main="p53/left")
myHeatmap(gaertner.right[p53,], main="p53/right")
myHeatmap(gaertner.left[wnt,], main="Wnt/left")
myHeatmap(gaertner.right[wnt,], main="Wnt/right")
myHeatmap(gaertner.left[fn1,], main="FN1 targets/left")
myHeatmap(gaertner.right[fn1,], main="FN1 targets/right")
myHeatmap(gaertner.left[acyl,], main="Acyl remodelling/left")
myHeatmap(gaertner.right[acyl,], main="Acyl remodelling/right")
myHeatmap(gaertner.left[g1s,], main="G1S transition/left")
myHeatmap(gaertner.right[g1s,], main="G1S transition/right")
dev.off()
plotGene <- function(eset,genesymbol,...) {
  ind <- mmatch(genesymbol, fData(eset)$GeneSymbol)[[1]]
  group <- factor(eset$group)
  sub <- colMeans(rowscale(exprs(eset)[ind,,drop=FALSE]))
  boxplot(sub~group, col=c("yellow", "red"), ylab="z-score", ...)
          
}
pdf(outfig("selected-genes.pdf"), width=4L, height=8L)
par(mfrow=c(2,1), mar=c(4,4,2,1), mgp=c(2,1,0))
plotGene(gaertner.left, "FOXO1", main="FOXO1/left")
plotGene(gaertner.right, "FOXO1", main="FOXO1/right")
plotGene(gaertner.left, "FOXO3", main="FOXO3/left")
plotGene(gaertner.right, "FOXO3", main="FOXO3/right")
plotGene(gaertner.left, "GSK3B", main="GSK3B/left")
plotGene(gaertner.right, "GSK3B", main="GSK3B/right")
dev.off()
