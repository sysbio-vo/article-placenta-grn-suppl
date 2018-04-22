library(cowplot)
library(ggfortify)
library(beadarray)
library(arrayQualityMetrics)
library(sva)
library(dplyr)
source("plots_utils.R")

### General variables

studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

## Choose between cohorts
i = which(studies$ID=="london")
#i = which(studies$ID=="oslo")

### Read data
path = paste("../raws/", studies[i,]$ID, sep="")
pdata = read.table(paste("../pdata/", studies[i,]$ID, "_pdata_untracked.tsv", sep=""), 
                   sep="\t", head=TRUE, stringsAsFactors = FALSE)
raw.data <- readIllumina(dir=path, sampleSheet=paste("../pdata/", studies[i,]$ID, "_pdata_untracked.csv", sep=""),
                     illuminaAnnotation="Humanv3")

### Normalize and QC

## P95 and P05 scores ratio
sID <- as.factor(pdata$Sentrix_ID)
ht12metrics <- c()
for (id in levels(sID)) {
  ht12metrics <- rbind(ht12metrics, read.table(paste(path, "/", id, "/Metrics.txt", sep=""),
                            sep = "\t", header = TRUE, as.is = TRUE))
}

pl <- illuPRatio(ht12metrics)
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PRatio.pdf", sep=""),
          base_height=3, base_aspect_ratio = 2.5, pl, ncol=1)
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PRatio.svg", sep=""),
          base_height=3, base_aspect_ratio = 2.5, pl, ncol=1)

## Array image plot
#imageplot(raw.data, array=16,high="white",low="darkgreen",zlim=c(4,10))
#combinedControlPlot(raw.data,array=1)

## Summarization and mormalization 
myMean = function(x) mean(x,na.rm=TRUE)
mySd = function(x) sd(x,na.rm=TRUE)
greenChannel = new("illuminaChannel", greenChannelTransform, illuminaOutlierMethod, myMean, mySd,"G")

eset.sum <- beadarray::summarize(raw.data, channelList = list(greenChannel))
boxplot(exprs(eset.sum),outline=FALSE)
eset.norm <- normaliseIllumina(eset.sum, method="neqc", transform="none")
boxplot(exprs(eset.norm),outline=FALSE)

exprs <- exprs(eset.norm)
if (studies[i,]$ID=="london") {
  third = which(pdata$Trimester!="Third")
} else {
  third = c(1:nrow(pdata))
}

#save(list=c("exprs", "pdata"), file="out")
#load("out")
eset = ExpressionSet(assayData=exprs[,third],
                     phenoData = AnnotatedDataFrame(pdata[third,]))

arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/QC/", studies[i,]$ID,"/AQM_report_", studies[i,]$ID, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))

exprs[which(is.na(exprs))] <- min(exprs[!is.na(exprs)])
#which(is.na(exprs), arr.ind = T)
pca = prcomp(t(exprs[,third]))
pl <- pcaPlots(pca, pdata[third,], c("Condition", "Sentrix_ID"))

## Save plot for manual quality control
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA.svg", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

## Remove outliers
exl = which(pdata$QC=="Outlier")
pdata <- pdata[-exl,]
exprs <- exprs[,-exl]
if (studies[i,]$ID=="london") {
  third = which(pdata$Trimester!="Third")
} else {
  third = c(1:nrow(pdata))
}

eset = ExpressionSet(assayData=exprs[,third],
                     phenoData = AnnotatedDataFrame(pdata[third,]))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/QC/", studies[i,]$ID,"/AQM_report_", studies[i,]$ID, "_nooutliers", sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))


### Batch-effect removal

if (studies[i,]$ID=="london") {
  batch = as.factor(pdata$Batch)
} else {
  batch = as.factor(pdata$Sentrix_ID)
}
mod = model.matrix(~as.factor(Condition), data=pdata)
exprs.nobatch = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

## Perform PCA and create plots
pca = prcomp(t(exprs[,third]))
pca.nobatch = prcomp(t(exprs.nobatch[,third]))
pl <- pcaPlots(pca, pdata[third,], c("Condition", "Sentrix_ID"))
pl.nobatch <- pcaPlots(pca.nobatch, pdata[third,], c("Condition", "Sentrix_ID"))

## Save plot for manual quality control
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers.svg", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]])
save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch.svg", sep=""),
          base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]])

if (studies[i,]$ID=="oslo") {
  pca.nobatch = prcomp(t(exprs.nobatch[,third]))
  pl.nobatch <- pcaPlots(pca.nobatch, pdata[third,], c("Condition", "ConditionDetailed", "Onset", "acogPE"), ncol=2)
  save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch_cond.pdf", sep=""),
            base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]], nrow=2)
  save_plot(paste("../plots/QC/", studies[i,]$ID, "/", studies[i,]$ID, "_PCA_nooutliers_nobatch_cond.svg", sep=""),
            base_height=3, base_aspect_ratio = pl.nobatch[[2]], pl.nobatch[[1]], nrow=2)
}

eset = ExpressionSet(assayData=exprs.nobatch[,third],
                     phenoData = AnnotatedDataFrame(pdata[third,]))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/QC/", studies[i,]$ID,"/AQM_report_", studies[i,]$ID, "_nooutliers_nobatch", sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))

### Save all the data (without outliers)
write.table(exprs.nobatch, paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""), sep="\t", quote=FALSE)
