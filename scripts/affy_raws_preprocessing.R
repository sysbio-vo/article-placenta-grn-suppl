# Load hgu133 BrainArray packages
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
library(affy)
library(cowplot)
library(ggfortify)
library(arrayQualityMetrics)
library(stringr) 
library(sva)
source("plots_utils.R")

### General variables

studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

## Choose between cohorts
i = which(studies$ID=="pittsburgh")

path = paste("../raws/", studies[i,]$ID, sep="")
pdata = read.table(paste("../pdata/", studies[i,]$ID, "_pdata_untracked.tsv", sep=""), 
                   sep="\t", head=TRUE, stringsAsFactors = FALSE)
pd <- AnnotatedDataFrame(pdata)

### Read and normalize the data

affyData = ReadAffy(phenoData=pd, sampleNames=pd$SampleAccessionNumber, filenames=as.character(rownames(pd)),
                    celfile.path=path)

affyData@cdfName <- "hgu133plus2hsentrezgcdf"
eset = rma(affyData)

### Processing

exprs <- exprs(eset)

## Perform PCA and create plots
pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pdata, c("Condition"))

## Save plot for manual quality control
save_plot(paste("../plots/QC/", studies[i,]$ID, "_PCA.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]]/2, pl[[1]])
save_plot(paste("../plots/QC/", studies[i,]$ID, "_PCA.svg", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]]/2, pl[[1]])

### Manual QC step

## Check the plots and decide if the batch correction is necessary
## Proceed with outlier detection

### arrayQualityMetrics report

## We use only one dataset, original or nobatch, decided during manual step
arrayQualityMetrics(expressionset = eset,
                   outdir = paste("../plots/QC/AQM_report_", studies[i,]$ID, sep=""),
                   force = TRUE,
                   do.logtransform = FALSE,
                   intgroup = c("Condition"))

arrayQualityMetrics(expressionset = affyData,
                    outdir = paste("../plots/QC/AQM_report_raw_", studies[i,]$ID, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))

## Save all the data, e.g. with outliers, just in case, the main idea is to add the info to pdata during manual steps
write.table(exprs, paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""), sep="\t", quote=FALSE)
