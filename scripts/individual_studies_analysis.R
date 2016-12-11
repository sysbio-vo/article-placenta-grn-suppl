library(illuminaHumanv3.db)
library(WGCNA)
library(limma)

### General variables

studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

## Choose between cohorts
i = which(studies$ID=="london")
#i = which(studies$ID=="oslo")

pdata <- read.table(paste("../pdata/", studies[i,]$ID, "_pdata.tsv", sep=""), 
                   sep="\t", head=TRUE, stringsAsFactors = FALSE)
exprs <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""),
                   header=TRUE, check.names=FALSE)

## Exclude outliers
exl = which(pdata$QC=="Outlier")
pdata <- pdata[-exl,]
exprs <- exprs[,-exl]

## In London cohort we have samples from Oslo, exclude them for now
if (studies[i,]$ID=="london") {
  exl = which(pdata$Trimester==3)
  pdata <- pdata[-exl,]
  exprs <- exprs[,-exl]
}

### Getting probe to gene unique correspondence

## Get probeset to entrezid mapping
probesetsID <- rownames(exprs)
probesetsID_EntrezID<-select(get(paste(studies[i,]$platformAbbr, ".db", sep="")), probesetsID, "ENTREZID")

## Replace probesetsIDs with gene IDs in expression data matrix

# Exclude NA probesets
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
# Exclude probesets mapped to different genes simultaneously
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
# Filter expression matrix based on left probesets
exprs <- exprs[which(rownames(exprs) %in% probesetsID_EntrezID$PROBEID),]

# Select one probeset among the probesets mapped to the same gene based on maximum average value across the samples
collapsed = collapseRows(exprs, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
exprs.unique <-collapsed$datETcollapsed

### Differentially expressed genes

## Find degs using linear models
design = model.matrix(~as.factor(Condition), data=pdata)
colnames(design) <- c("Intercept", "LOWvsHIGH")
fit <- lmFit(exprs.unique, design)
fit <- eBayes(fit)
degs <- topTable(fit, coef="LOWvsHIGH", adjust.method="fdr", number=nrow(fit))

## Filter by adj.P.Val
degs <- degs[degs$adj.P.Val<0.05,]


