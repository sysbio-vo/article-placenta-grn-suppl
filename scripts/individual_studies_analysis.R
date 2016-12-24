library(illuminaHumanv3.db)
library(WGCNA)
library(limma)
source("degs_utils.R")

### General variables
TEST = FALSE
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

exprs.unique <- getUniqueProbesets(exprs, studies[i,]$platformAbbr)

### Differentially expressed genes

if (studies[i,]$ID=="london") {
  degs <- getDEGS(c("Low risk", "High risk"), pdata, exprs.unique)
} else {
  degs <- getDEGS(c("Control", "Preeclampsia"), pdata, exprs.unique)
}

write.table(degs, paste("../degs/", studies[i,]$ID, "_degs.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)

# Filter by p-value and logFC, write to file short list
if (studies[i,]$ID=="london") {
  degs <- filterDEGS(degs, 0.05, 0.7, adj=FALSE)
} else {
  degs <- filterDEGS(degs, 0.05, 0.7)
}

write.table(degs, paste("../degs/", studies[i,]$ID, "_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
