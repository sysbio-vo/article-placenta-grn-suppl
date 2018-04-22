library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
library(affy)
library(illuminaHumanv3.db)
library(WGCNA)
library(limma)
library(sva)
library(beadarray)
library(arrayQualityMetrics)
source("plots_utils.R")
source("degs_utils.R")

### General variables
studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

## Read from files
i = which(studies$ID=="london")
pdata.london <- read.table(paste("../pdata/", studies[i,]$ID, "_pdata_untracked.tsv", sep=""), 
                    sep="\t", head=TRUE, stringsAsFactors = FALSE)
exprs.london <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""),
                    header=TRUE, check.names=FALSE)
i = which(studies$ID=="oslo")
pdata.oslo <- read.table(paste("../pdata/", studies[i,]$ID, "_pdata_untracked.tsv", sep=""), 
                           sep="\t", head=TRUE, stringsAsFactors = FALSE)
exprs.oslo <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""),
                           header=TRUE, check.names=FALSE)
i = which(studies$ID=="pittsburgh")
pdata.pitt <- read.table(paste("../pdata/", studies[i,]$ID, "_pdata_untracked.tsv", sep=""), 
                         sep="\t", head=TRUE, stringsAsFactors = FALSE)
exprs.pitt <- read.table(paste("../exprs/", studies[i,]$ID, "_exprs.tsv", sep=""),
                         header=TRUE, check.names=FALSE)

## Merge Oslo and London into one dataset
cols <- c("SampleID", "Condition", "Trimester", "StudyID", "QC", "ReplicateSample")
pdata <- rbind(pdata.london[,cols], pdata.oslo[,cols])
exprs <- merge(exprs.london, exprs.oslo, by="row.names")
rownames(exprs) <- exprs[,1]
exprs <- exprs[,-1]

## Exclude outliers
exl = which(pdata$QC=="Outlier")
pdata <- pdata[-exl,]

# PCA before cross-experiment batch removal

pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "StudyID", "ReplicateSample"), ncol=2)
save_plot("../plots/Merge/LO_integrated_PCA.pdf",
          base_height=4, base_aspect_ratio = pl[[2]]/2, pl[[1]], nrow=2)
save_plot("../plots/Merge/LO_integrated_PCA.svg",
          base_height=4, base_aspect_ratio = pl[[2]]/2, pl[[1]], nrow=2)

### Remove batch caused by two different studies

batch = as.factor(pdata$StudyID)
mod = model.matrix(~as.factor(Trimester), data=pdata)
exprs.nobatch = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

## PCA after batch-effect removal
pca = prcomp(t(exprs.nobatch))
pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "StudyID", "ReplicateSample"), ncol=2)
save_plot("../plots/Merge/LO_integrated_PCA_nobatch.pdf",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.5, pl[[1]], nrow=2)
save_plot("../plots/Merge/LO_integrated_PCA_nobatch.svg",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.5, pl[[1]], nrow=2)


## Separate datasets
# London
l <- rownames(pdata[which(pdata$StudyID=="London" & pdata$Trimester=="First"),])
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata[which(rownames(pdata) %in% l),],
               c("Condition", "Trimester"))
save_plot("../plots/Merge/LO_london_integrated_PCA_nobatch.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LO_london_integrated_PCA_nobatch.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

# Oslo
l <- rownames(pdata[which(pdata$Trimester=="Third"),])
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata[which(rownames(pdata) %in% l),],
               c("Condition", "ReplicateSample"))
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

### Getting probe to gene unique correspondence

exprs.unique <- getUniqueProbesets(exprs.nobatch, studies[1,]$platformAbbr)

## Check if probesets elimintation distorted PCA plot (didn't change much)
pca = prcomp(t(exprs.unique))
pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "StudyID", "ReplicateSample"), ncol=2)
save_plot("../plots/Merge/LO_integrated_PCA_nobatch_unique_probesets.pdf",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.5, pl[[1]], nrow=2)
save_plot("../plots/Merge/LO_integrated_PCA_nobatch_unique_probesets.svg",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.5, pl[[1]], nrow=2)

## Separate datasets
# London
l <- rownames(pdata[which(pdata$StudyID=="London" & pdata$Trimester=="First"),])
pca = prcomp(t(exprs.unique[,which(colnames(exprs.unique) %in% l)]))
pl <- pcaPlots(pca, pdata[which(rownames(pdata) %in% l),],
               c("Condition", "Trimester"))
save_plot("../plots/Merge/LO_london_integrated_PCA_nobatch_unique_probesets.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LO_london_integrated_PCA_nobatch_unique_probesets.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

# Oslo
l <- rownames(pdata[which(pdata$Trimester=="Third"),])
pca = prcomp(t(exprs.unique[,which(colnames(exprs.unique) %in% l)]))
pl <- pcaPlots(pca, pdata[which(rownames(pdata) %in% l),],
               c("Condition", "ReplicateSample"))
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch_unique_probesets.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch_unique_probesets.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

## Merge Oslo+London and Pittsburgh into one dataset

cols <- c("SampleID", "Condition", "Trimester", "StudyID", "Batch")
pdata$Batch <- 1
pdata.pitt$Batch <- 2
pdata.all <- rbind(pdata[,cols], pdata.pitt[,cols])
pdata.all$ConditionDetailed <- pdata.all$Condition
pdata.all$Condition[which(pdata.all$Condition=="Normal RI")] <- "Norma"

# Get EntrezID from Affy probe IDs
probesetsID_EntrezID<-select(hgu133plus2hsentrezg.db, rownames(exprs.pitt), "ENTREZID")
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]

exprs.pitt <- exprs.pitt[which(rownames(exprs.pitt) %in% probesetsID_EntrezID$PROBEID),]
rownames(exprs.pitt) <- probesetsID_EntrezID[match(rownames(exprs.pitt), probesetsID_EntrezID$PROBEID),]$ENTREZID

exprs <- merge(exprs.unique, exprs.pitt, by="row.names")
rownames(exprs) <- exprs[,1]
exprs <- exprs[,-1]


# PCA before cross-experiment batch removal

pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pdata.all, c("Condition", "Trimester", "StudyID", "Batch"), ncol=2)
save_plot("../plots/Merge/LOP_integrated_PCA.pdf",
          base_height=4, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2)
save_plot("../plots/Merge/LOP_integrated_PCA.svg",
          base_height=4, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2)

### Remove batch caused by two different studies

batch = as.factor(pdata.all$Batch)
#pdata.all$mod <- paste(pdata.all$Trimester, pdata.all$Condition)
mod = model.matrix(~as.factor(Trimester), data=pdata.all)
exprs.nobatch = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

## PCA after batch-effect removal
pca = prcomp(t(exprs.nobatch))
pl <- pcaPlots(pca, pdata.all, c("Condition", "Trimester", "StudyID", "Batch"), ncol=2)
save_plot("../plots/Merge/LOP_integrated_PCA_nobatch.pdf",
          base_height=3.3, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2)
save_plot("../plots/Merge/LOP_integrated_PCA_nobatch.svg",
          base_height=3.3, base_aspect_ratio = pl[[2]], pl[[1]], nrow=2)

## Separate datasets
# Trimester
l <- pdata.all[which(pdata.all$Trimester=="First"),]$SampleID
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))

pl <- pcaPlots(pca, pdata.all[which(pdata.all$SampleID %in% l),],
               c("Condition", "ConditionDetailed", "StudyID"), ncol=3)
save_plot("../plots/Merge/LOP_first_integrated_PCA.pdf",
          base_height=3, base_aspect_ratio = pl[[2]]*1.2, pl[[1]])
save_plot("../plots/Merge/LOP_first_integrated_PCA.svg",
          base_height=3, base_aspect_ratio = pl[[2]]*1.2, pl[[1]])

l <- pdata.all[which(pdata.all$Trimester=="Third"),]$SampleID
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata.all[which(pdata.all$SampleID %in% l),],
               c("Condition", "StudyID"))
save_plot("../plots/Merge/LOP_third_integrated_PCA.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LOP_third_integrated_PCA.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

# London
l <- rownames(pdata.all[which(pdata.all$StudyID=="London" & pdata.all$Trimester=="First"),])
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata.all[which(rownames(pdata.all) %in% l),],
               c("Condition", "Trimester"))
save_plot("../plots/Merge/LOP_london_integrated_PCA.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LOP_london_integrated_PCA.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])


write.table(exprs.nobatch, paste("../exprs/exprs_all.tsv", sep=""), sep="\t", quote=FALSE)
rownames(pdata.all) <- pdata.all$SampleID
col <- c("SampleID", "Condition", "ConditionDetailed", "Trimester", "StudyID")
write.table(pdata.all[, col], paste("../pdata/pdata_all.tsv", sep=""), sep="\t", quote=FALSE)
