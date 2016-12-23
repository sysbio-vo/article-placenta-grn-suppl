library(illuminaHumanv3.db)
library(WGCNA)
library(limma)
library(sva)
library(arrayQualityMetrics)
source("pcaPlots.R")
library(beadarray)

### General variables

studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

## Read from files
i = which(studies$ID=="london")
pdata.london <- read.table(paste("../pdata/", studies[i,]$ID, "_pdata.tsv", sep=""), 
                    sep="\t", head=TRUE, stringsAsFactors = FALSE)
exprs.london <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""),
                    header=TRUE, check.names=FALSE)
i = which(studies$ID=="oslo")
pdata.oslo <- read.table(paste("../pdata/", studies[i,]$ID, "_pdata.tsv", sep=""), 
                           sep="\t", head=TRUE, stringsAsFactors = FALSE)
exprs.oslo <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""),
                           header=TRUE, check.names=FALSE)

## Merge into one dataset
pdata <- rbind(pdata.london, pdata.oslo)
exprs <- merge(exprs.london, exprs.oslo, by="row.names")
rownames(exprs) <- exprs[,1]
exprs <- exprs[,-1]

## Exclude outliers
exl = which(pdata$QC=="Outlier")
pdata <- pdata[-exl,]
exprs <- exprs[,-exl]

### Remove batch caused by two different studies

batch = as.factor(pdata$Study_ID)
mod = model.matrix(~as.factor(Condition), data=pdata)
combat_edata = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs.nobatch <- combat_edata

## Check with PCA plot

## Integrated dataset
pca = prcomp(t(exprs.nobatch))
pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "Study_ID", "Sample_Alt_Name"), "Integrated data", ncol=2)
save_plot("../plots/PCA/integrated_PCA_nobatch.pdf",
          base_height=5.5, base_aspect_ratio = 1.6, pl)

eset = ExpressionSet(assayData=as.matrix(exprs.nobatch), phenoData = AnnotatedDataFrame(pdata))
arrayQualityMetrics(expressionset = eset,
                    outdir = "../plots/AQM/AQM_report_nobatch_integrated",
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))


## Separate datasets
# London
l <- rownames(pdata[which(pdata$Study_ID=="london"),])
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata[which(pdata$Study_ID=="london"),],
               c("Condition", "Trimester"), "Integrated data. London")
save_plot("../plots/PCA/london_integrated_PCA_nobatch.pdf",
          base_height=4, base_aspect_ratio = 2.8, pl)

# Oslo
l <- rownames(pdata[which(pdata$Study_ID=="oslo"),])
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata[which(pdata$Study_ID=="oslo"),],
               c("Condition", "Trimester"), "Integrated data. London")
save_plot("../plots/PCA/oslo_integrated_PCA_nobatch.pdf",
          base_height=4, base_aspect_ratio = 2.8, pl)

### Getting probe to gene unique correspondence

exprs <- exprs.nobatch
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

## Check if probesets elimintation distorted PCA plot (didn't change much)
pca = prcomp(t(exprs.unique))
pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "Study_ID", "Sample_Alt_Name"), "Integrated data", ncol=2)
save_plot("../plots/PCA/integrated_PCA_nobatch_unique_probesets.pdf",
          base_height=5.5, base_aspect_ratio = 1.6, pl)

### Differentially expressed genes

## Find degs using linear models

# High vs Preeclampsia
#design = model.matrix(~as.factor(Condition), data=pdata)

# Low vs Preeclampsia

# High vs Control
