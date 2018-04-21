library(beadarray)
library(limma)
source("plots_utils.R")

studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

## Choose between cohorts
i = which(studies$ID=="london")
#i = which(studies$ID=="oslo")

pdata = read.table(paste("../pdata/", studies[i,]$ID, "_pdata_untracked.tsv", sep=""), 
                   sep="\t", head=TRUE, stringsAsFactors = FALSE)

write.table(pdata, file=paste("../pdata/", studies[i,]$ID, "_pdata_untracked.csv", sep=""),
            sep=",", quote=FALSE, row.names=FALSE)


path = paste("../raws/", studies[i,]$ID, sep="")
data <- readIllumina(dir=path,sampleSheet=paste("../pdata/", studies[i,]$ID, "_pdata_untracked.csv", sep=""),
                     illuminaAnnotation="Humanv3")

dev.off()
imageplot(data, array=1,high="darkgreen",low="lightgreen",zlim=c(4,10))
imageplot(data, array=2,high="darkgreen",low="lightgreen",zlim=c(4,10))
imageplot(data, array=3,high="darkgreen",low="lightgreen",zlim=c(4,10))

combinedControlPlot(data,array=1)
combinedControlPlot(data,array=2)
combinedControlPlot(data,array=3)

eset.ill <- summarize(data)

boxplot(exprs(eset.ill),outline=FALSE)

eset.norm <- normaliseIllumina(eset.ill)

boxplot(exprs(eset.norm),outline=FALSE)

exprs <- exprs(eset.norm)
exl = which(pdata$Trimester=="Third")
pdata <- pdata[-exl,]
exprs <- exprs[,-exl]

eset = ExpressionSet(assayData=exprs, phenoData = AnnotatedDataFrame(pdata))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/QC/AQM_report_alt_first_", studies[i,]$ID, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Sample_Group"))

exl = which(pdata$QC=="Outlier")
exl = 16
pdata <- pdata[-exl,]
exprs <- exprs[,-exl]
exprs.nobatch <- exprs.nobatch[,-exl]
exprs.unique <- exprs.unique[,-exl]

eset = ExpressionSet(assayData=exprs, phenoData = AnnotatedDataFrame(pdata))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/QC/AQM_report_alt_nooutliers_", studies[i,]$ID, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Sample_Group"))

batch = as.factor(pdata$Batch)
mod = model.matrix(~as.factor(Sample_Group), data=pdata)
exprs.nobatch = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

exprs.nobatch[which(is.na(exprs.nobatch))] <- 0
pca.nobatch = prcomp(t(exprs.nobatch))
pl.nobatch <- pcaPlots(pca.nobatch, pdata, c("Sample_Group", "Sentrix_ID"))

## Save plot for manual quality control
save_plot(paste("../plots/QC/", studies[i,]$ID, "_alt_PCA.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl.nobatch[[1]])
save_plot(paste("../plots/QC/", studies[i,]$ID, "_alt_PCA.svg", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl.nobatch[[1]])

exprs.unique <- getUniqueProbesets(exprs.nobatch, studies[1,]$platformAbbr)
pca.unique = prcomp(t(exprs.unique))
pl.unique <- pcaPlots(pca.unique, pdata, c("Sample_Group", "Sentrix_ID"))

## Save plot for manual quality control
save_plot(paste("../plots/QC/", studies[i,]$ID, "_altu_PCA.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl.unique[[1]])
save_plot(paste("../plots/QC/", studies[i,]$ID, "_altu_PCA.svg", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]], pl.unique[[1]])

pdata$Condition <- pdata$Sample_Group
degs <- getDEGS(c("Normal RI", "High RI"), pdata, exprs.unique)
degs <- filterDEGS(degs, 0.05, 0.5, adj=FALSE)

mmdegs <- degs
