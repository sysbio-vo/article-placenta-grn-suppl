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
l <- rownames(pdata[which(pdata$StudyID=="Oslo" & pdata$Trimester=="Third"),])
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata[which(rownames(pdata) %in% l),],
               c("Condition", "Trimester"))
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

### Getting probe to gene unique correspondence

exprs.unique <- getUniqueProbesets(exprs.nobatch, studies[1,]$platformAbbr)
#write.table(exprs.unique, "../exprs/exprs_all.tsv", sep="\t", row.names = TRUE, quote=FALSE)


## Check if probesets elimintation distorted PCA plot (didn't change much)
pca = prcomp(t(exprs.unique))
pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "StudyID", "ReplicateSample"), ncol=2)
save_plot("../plots/Merge/LO_integrated_PCA_nobatch_unique_probesets.pdf",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.1, pl[[1]], nrow=2)
save_plot("../plots/Merge/LO_integrated_PCA_nobatch_unique_probesets.svg",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.1, pl[[1]], nrow=2)

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
l <- rownames(pdata[which(pdata$StudyID=="Oslo" & pdata$Trimester=="Third"),])
pca = prcomp(t(exprs.unique[,which(colnames(exprs.unique) %in% l)]))
pl <- pcaPlots(pca, pdata[which(rownames(pdata) %in% l),],
               c("Condition", "Trimester"))
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch_unique_probesets.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/Merge/LO_oslo_integrated_PCA_nobatch_unique_probesets.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
