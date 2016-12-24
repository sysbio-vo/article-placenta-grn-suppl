library(illuminaHumanv3.db)
library(WGCNA)
library(limma)
library(sva)
library(arrayQualityMetrics)
source("pcaPlots.R")
source("degs_utils.R")
library(beadarray)

### General variables
TEST = FALSE
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
if (TEST) {
  
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

}

### Getting probe to gene unique correspondence

exprs.unique <- getUniqueProbesets(exprs.nobatch, studies[i,]$platformAbbr)

## Check if probesets elimintation distorted PCA plot (didn't change much)
if (TEST) {
  pca = prcomp(t(exprs.unique))
  pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "Study_ID", "Sample_Alt_Name"), "Integrated data", ncol=2)
  save_plot("../plots/PCA/integrated_PCA_nobatch_unique_probesets.pdf",
          base_height=5.5, base_aspect_ratio = 1.6, pl)
}

### Differentially expressed genes

## Preeclampsia versus high and low

# High vs Preeclampsia
high.degs <- getDEGS(c("High risk", "Preeclampsia"), pdata, exprs.unique)
write.table(high.degs, paste("../degs/HP_degs.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
high.degs <- filterDEGS(high.degs, 0.01, 3.5)

# Low vs Preeclampsia
low.degs <- getDEGS(c("Low risk", "Preeclampsia"), pdata, exprs.unique)
write.table(low.degs, paste("../degs/LP_degs.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
low.degs <- filterDEGS(low.degs, 0.01, 3.5)

# Some comparison
common <- intersect(low.degs$SYMBOL, high.degs$SYMBOL)
Unique_Gene <- !(low.degs$SYMBOL %in% common)
low.degs <- cbind(Unique_Gene, low.degs)
Unique_Gene <- !(high.degs$SYMBOL %in% common)
high.degs <- cbind(Unique_Gene, high.degs)

write.table(high.degs, paste("../degs/HP_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
write.table(low.degs, paste("../degs/LP_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)

## Control versus high and low

# High vs Control

high.degs <- getDEGS(c("High risk", "Control"), pdata, exprs.unique)
write.table(high.degs, paste("../degs/HC_degs.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
high.degs <- filterDEGS(high.degs, 0.01, 3.5)

# Low vs Control

low.degs <- getDEGS(c("Low risk", "Control"), pdata, exprs.unique)
write.table(low.degs, paste("../degs/LC_degs.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
low.degs <- filterDEGS(low.degs, 0.01, 3.5)

# Some comparison

common <- intersect(low.degs$SYMBOL, high.degs$SYMBOL)
Unique_Gene <- !(low.degs$SYMBOL %in% common)
low.degs <- cbind(Unique_Gene, low.degs)
Unique_Gene <- !(high.degs$SYMBOL %in% common)
high.degs <- cbind(Unique_Gene, high.degs)

write.table(high.degs, paste("../degs/HC_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
write.table(low.degs, paste("../degs/LC_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)

## Control vs Preeclampsia

degs <- getDEGS(c("Control", "Preeclampsia"), pdata, exprs.unique)
write.table(degs, paste("../degs/CP_degs.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
degs <- filterDEGS(degs, 0.05, 0.7)
write.table(degs, paste("../degs/CP_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)

# Some comparison
degs.int <- read.table(paste("../degs/CP_degs_short.tsv", sep=""),
                           header=TRUE, check.names=FALSE, sep = "\t")

degs.ind <- read.table(paste("../degs/oslo_degs_short.tsv", sep=""),
                       header=TRUE, check.names=FALSE, sep = "\t")

common <- intersect(degs.int$SYMBOL, degs.ind$SYMBOL)
Unique_Gene <- !(degs.int$SYMBOL %in% common)
degs.int <- cbind(Unique_Gene, degs.int)
Unique_Gene <- !(degs.ind$SYMBOL %in% common)
degs.ind <- cbind(Unique_Gene, degs.ind)

write.table(degs.int, paste("../degs/CP_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
write.table(degs.ind, paste("../degs/oslo_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)
