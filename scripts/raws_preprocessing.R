library(cowplot)
library(ggfortify)
library(beadarray)
library(arrayQualityMetrics)
library(sva)
source("pcaPlots.R")

### General variables

studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

## Choose between cohorts
i = which(studies$ID=="london")
#i = which(studies$ID=="oslo")

path = paste("../raws/", studies[i,]$ID, sep="")
pdata = read.table(paste("../pdata/", studies[i,]$ID, "_pdata.tsv", sep=""), 
                   sep="\t", head=TRUE, stringsAsFactors = FALSE)
idatfiles = dir(path, pattern="idat", full.name=TRUE)

### Read and normalize the data

raw.data <- readIdatFiles(idatfiles)
raw.data@phenoData = AnnotatedDataFrame(pdata)
#N.data <- normaliseIllumina(channel(raw.data, "Green"), method="neqc", transform="none")
N.data <- normaliseIllumina(channel(raw.data, "Green"), method="quantile", transform="log2")
## In case of quantile do the following
exprs <- exprs(N.data)
# Remove all rows with non-finite values
exprs <- exprs[!rowSums(!is.finite(exprs)),]
# Replace all non-finite values with 0
exprs[!is.finite(exprs)] <- 0
# Change back
N.data@assayData$exprs <- exprs

### Processing

## Batch-effect removal
batch = as.factor(pData(N.data)$Sentrix_ID)
mod = model.matrix(~as.factor(Condition), data=pData(N.data))
combat_edata = ComBat(dat=exprs(N.data), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
N.data.nobatch <- N.data
N.data.nobatch@assayData$exprs <- combat_edata

## In case of quantile do the following
exprs <- exprs(N.data.nobatch)
# Remove all rows with non-finite values
exprs <- exprs[!rowSums(!is.finite(exprs)),]
# Replace all non-finite values with 0
exprs[!is.finite(exprs)] <- 0
# Change back
N.data.nobatch@assayData$exprs <- exprs

## Extract expression data
N.exprs <- exprs(N.data)
N.exprs.nobatch <- exprs(N.data.nobatch)

## In London cohort we have samples from Oslo, exclude them to plot PCA properly
if (studies[i,]$ID=="london") {
  exl = which(pdata$Trimester==3)
  pdata <- pdata[-exl,]
  N.exprs <- N.exprs[,-exl]
  N.exprs.nobatch <- N.exprs.nobatch[,-exl]
}

## Perform PCA and create plots
pca = prcomp(t(N.exprs))
pca.nobatch = prcomp(t(N.exprs.nobatch))
pl <- pcaPlots(pca, pdata, c("Condition", "Sentrix_ID"), studies[i,]$ID)
pl.nobatch <- pcaPlots(pca.nobatch, pdata, c("Condition", "Sentrix_ID"), studies[i,]$ID)

## Save plot for manual quality control
save_plot(paste("../plots/PCA/", studies[i,]$ID, "_preliminary_PCA.pdf", sep=""),
          base_height=4, base_aspect_ratio = 2.6, pl)
save_plot(paste("../plots/PCA/", studies[i,]$ID, "_preliminary_PCA_nobatch.pdf", sep=""),
          base_height=4, base_aspect_ratio = 2.6, pl.nobatch)

### Manual QC step

## Check the plots and decide if the batch correction is necessary
## Proceed with outlier detection

### arrayQualityMetrics report

## We use only one dataset, original or nobatch, decided during manual step
eset = ExpressionSet(assayData=N.exprs.nobatch, phenoData = AnnotatedDataFrame(pdata))
arrayQualityMetrics(expressionset = eset,
                   outdir = paste("../plots/AQM/AQM_report_nobatch_", studies[i,]$ID, sep=""),
                   force = TRUE,
                   do.logtransform = FALSE,
                   intgroup = c("Condition"))

### Manual outlier detection step
## Read AQM report and decide, which samples should be excluded, edit pdata file

### arrayQualityMetrics report without outliers
exl = which(pdata$QC=="Outlier")
pdata <- pdata[-exl,]
N.exprs.nobatch <- N.exprs.nobatch[,-exl]
eset = ExpressionSet(assayData=N.exprs.nobatch, phenoData = AnnotatedDataFrame(pdata))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/AQM/AQM_report_nobatch_nooutliers_", studies[i,]$ID, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))
## Final PCA plot
pca.final = prcomp(t(N.exprs.nobatch))
pl.final <- pcaPlots(pca.final, pdata, c("Condition", "Sentrix_ID"), studies[i,]$ID)
save_plot(paste("../plots/PCA/", studies[i,]$ID, "_final_PCA.pdf", sep=""),
          base_height=4, base_aspect_ratio = 2.6, pl.final)

### If the result is satisfying save the expression data for futher analysis

## Save all the data, e.g. with outliers, just in case, the main idea is to add the info to pdata during manual steps
N.exprs.nobatch <- exprs(N.data.nobatch)
write.table(N.exprs.nobatch, paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""), sep="\t", quote=FALSE)
