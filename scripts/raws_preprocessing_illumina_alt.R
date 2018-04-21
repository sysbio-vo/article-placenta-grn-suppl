library(beadarray)
library(sva)
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

i = which(studies$ID=="pittsburgh")
pdata.pitt <- read.table(paste("../pdata/", studies[i,]$ID, "_pdata_untracked.tsv", sep=""), 
                         sep="\t", head=TRUE, stringsAsFactors = FALSE)
exprs.pitt <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""),
                         header=TRUE, check.names=FALSE)

probesetsID_EntrezID<-select(hgu133plus2hsentrezg.db, rownames(exprs.pitt), "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]

exprs.pitt <- exprs.pitt[which(rownames(exprs.pitt) %in% probesetsID_EntrezID$PROBEID),]
rownames(exprs.pitt) <- probesetsID_EntrezID[match(rownames(exprs.pitt), probesetsID_EntrezID$PROBEID),]$ENTREZID


exprs.pitt <- exprs.pitt[which(rownames(exprs.pitt) %in% rownames(exprs.unique)),]
exprs.unique <- exprs.unique[which(rownames(exprs.unique) %in% rownames(exprs.pitt)),]
exprs.unique <- exprs.unique[order(match(rownames(exprs.unique), rownames(exprs.pitt))),]
exprs.all <- cbind(exprs.unique, exprs.pitt)

pd <- pdata[,c("ConditionSimple", "Sample_Name", "StudyID")]
colnames(pd) <- c("Condition", "SampleID", "StudyID")
pdpit <- pdata.pitt[, c("Condition", "SampleID", "StudyID")]

pd.all <- rbind(pd, pdpit)
exprs.all <- exprs.all[,pd.all$SampleID]

batch = as.factor(pd.all$StudyID)
mod = model.matrix(~as.factor(Condition), data=pd.all)
combat_edata = ComBat(dat=exprs.all, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs.nobatch <- combat_edata

pca = prcomp(t(exprs.nobatch))
pl <- autoplot(pca, data = pd.all, colour="StudyID") +
  scale_color_manual(values=c("red", "blue", "green", "black", "yellow", "pink", "brown"))

rownames(pd.all) <- pd.all$SampleID
degs <- getDEGS(c("Norma", "Preeclampsia"), pd.all, exprs.nobatch)

degs <- filterDEGS(degs, 0.05, 0.5)
