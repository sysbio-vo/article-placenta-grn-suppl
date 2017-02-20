library(illuminaHumanv3.db)
library(WGCNA)
library(limma)
library(sva)
library(arrayQualityMetrics)
source("plots_utils.R")
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

write.table(pdata, "../pdata/pdata.tsv", sep="\t", row.names = TRUE, quote=FALSE)

## Check PCA before batch correction
if (TEST) {
  pca = prcomp(t(exprs))
  pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "Study_ID", "Replicate_Sample"), ncol=2)
  save_plot("../plots/PCA/integrated_PCA.pdf",
            base_height=4, base_aspect_ratio = pl[[2]]/2, pl[[1]], nrow=2)
  save_plot("../plots/PCA/integrated_PCA.svg",
            base_height=4, base_aspect_ratio = pl[[2]]/2, pl[[1]], nrow=2)
}

### Remove batch caused by two different studies

batch = as.factor(pdata$Study_ID)
mod = model.matrix(~as.factor(Condition), data=pdata)
combat_edata = ComBat(dat=exprs, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
exprs.nobatch <- combat_edata

## Check with PCA plot
if (TEST) {
  
## Integrated dataset
pca = prcomp(t(exprs.nobatch))
pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "Study_ID", "Replicate_Sample"), ncol=2)
save_plot("../plots/PCA/integrated_PCA_nobatch.pdf",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.1, pl[[1]], nrow=2)
save_plot("../plots/PCA/integrated_PCA_nobatch.svg",
          base_height=3.3, base_aspect_ratio = pl[[2]]/1.1, pl[[1]], nrow=2)

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
               c("Condition", "Trimester"))
save_plot("../plots/PCA/london_integrated_PCA_nobatch.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/PCA/london_integrated_PCA_nobatch.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

# Oslo
l <- rownames(pdata[which(pdata$Study_ID=="oslo"),])
pca = prcomp(t(exprs.nobatch[,which(colnames(exprs.nobatch) %in% l)]))
pl <- pcaPlots(pca, pdata[which(pdata$Study_ID=="oslo"),],
               c("Condition", "Trimester"))
save_plot("../plots/PCA/oslo_integrated_PCA_nobatch.pdf",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])
save_plot("../plots/PCA/oslo_integrated_PCA_nobatch.svg",
          base_height=3, base_aspect_ratio = pl[[2]], pl[[1]])

}

### Getting probe to gene unique correspondence

exprs.unique <- getUniqueProbesets(exprs.nobatch, studies[i,]$platformAbbr)
if (TEST) {
  write.table(exprs.unique, "../exprs/exprs_all.tsv", sep="\t", row.names = TRUE, quote=FALSE)
}


## Check if probesets elimintation distorted PCA plot (didn't change much)
if (TEST) {
  pca = prcomp(t(exprs.unique))
  pl <- pcaPlots(pca, pdata, c("Condition", "Trimester", "Study_ID", "Replicate_Sample"), ncol=2)
  save_plot("../plots/PCA/integrated_PCA_nobatch_unique_probesets.pdf",
            base_height=3.3, base_aspect_ratio = pl[[2]]/1.1, pl[[1]], nrow=2)
  save_plot("../plots/PCA/integrated_PCA_nobatch_unique_probesets.svg",
            base_height=3.3, base_aspect_ratio = pl[[2]]/1.1, pl[[1]], nrow=2)
}

### Differentially expressed genes

# Low vs High
degs <- getDEGS(c("Low risk", "High risk"), pdata, exprs.unique)
write.table(degs, "../degs/LH_degs.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# High vs Preeclampsia
degs <- getDEGS(c("High risk", "Preeclampsia"), pdata, exprs.unique)
write.table(degs, "../degs/HP_degs.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# Low vs Preeclampsia
degs <- getDEGS(c("Low risk", "Preeclampsia"), pdata, exprs.unique)
write.table(degs, "../degs/LP_degs.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# High vs Control
degs <- getDEGS(c("High risk", "Control"), pdata, exprs.unique)
write.table(degs, "../degs/HC_degs.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# Low vs Control
degs <- getDEGS(c("Low risk", "Control"), pdata, exprs.unique)
write.table(degs, "../degs/LC_degs.tsv", sep="\t", row.names = FALSE, quote=FALSE)

## Control vs Preeclampsia
degs <- getDEGS(c("Control", "Preeclampsia"), pdata, exprs.unique)
write.table(degs, "../degs/CP_degs.tsv", sep="\t", row.names = FALSE, quote=FALSE)
degs <- filterDEGS(degs, 0.05, 0.7)
write.table(degs, paste("../degs/CP_degs_short.tsv", sep=""), sep="\t",
            row.names = FALSE, quote=FALSE)

## Comparison with old DEG

#CP
degs <- read.table("../degs/CP_degs.tsv", sep="\t", quote = "",
                   header=TRUE, check.names=FALSE)
degs <- filterDEGS(degs, 0.05, 0.7)
degs.old <- read.table(paste("../degs/oslo_old_degs_short.tsv", sep=""),
                        header=TRUE, check.names=FALSE, sep = "\t")
common <- intersect(degs$SYMBOL, degs.old$SYMBOL)

#LH
degs <- read.table("../degs/LH_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)
degs <- filterDEGS(degs, 0.01, 0.7, adj=FALSE)
degs.old <- read.table(paste("../degs/london_old_degs_short.tsv", sep=""),
                       header=TRUE, check.names=FALSE, sep = "\t")
degs.ind <- read.table(paste("../degs/london_degs_short.tsv", sep=""),
                       header=TRUE, check.names=FALSE, sep = "\t")
common <- intersect(degs$SYMBOL, degs.old$SYMBOL)
common <- intersect(degs$SYMBOL, degs.ind$SYMBOL)
