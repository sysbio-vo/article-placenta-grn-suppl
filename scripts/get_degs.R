library(illuminaHumanv3.db)
library(limma)
source("degs_utils.R")
source("plots_utils.R")

### General variables
studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

### Read data
exprs <- read.table("../exprs/exprs_all.tsv", sep="\t", quote = "",
                    header=TRUE, check.names=FALSE)
pdata  <- read.table("../pdata/pdata_all.tsv", sep="\t", quote = "",
                     header=TRUE, check.names=FALSE)
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

## Third trimester
ind <- which(pdata$Trimester=="Third")
degs <- getDEGS(c("Norma", "Preeclampsia"), pdata[ind,], exprs[,ind])
degs <- filterDEGS(degs, 0.05, 0.65)


exl = which(pdata.oslo$QC=="Outlier")
pdata.oslo <- pdata.oslo[-exl,]
exprs.unique <- getUniqueProbesets(exprs.oslo, studies[1,]$platformAbbr)
degs <- getDEGS(c("none", "early PE"), pdata.oslo, exprs.unique)
degs <- filterDEGS(degs, 0.05, 0.7)

## First trimester
ind <- which(pdata$Trimester=="First")
degs <- getDEGS(c("Norma", "Preeclampsia"), pdata[ind,], exprs[,ind])
degs1 <- filterDEGS(degs, 0.05, 0.6, adj=F)

ind <- which(pdata$Trimester=="First")
degs <- getDEGS(c("Norma", "High RI"), pdata[ind,], exprs[,ind])
degs2 <- filterDEGS(degs, 0.05, 0.6, adj=F)

ind <- which(pdata$Trimester=="First")
degs <- getDEGS(c("High RI", "Preeclampsia"), pdata[ind,], exprs[,ind])
degs2 <- filterDEGS(degs, 0.05, 0.6, adj=F)

exl = which(pdata.london$QC=="Outlier")
pdata.london <- pdata.london[-exl,]
exl = which(pdata.london$Trimester=="Third")
exprs.london <- exprs.london[,-exl]
exprs.unique <- getUniqueProbesets(exprs.london, studies[1,]$platformAbbr)
degs <- getDEGS(c("Normal RI", "High RI"), pdata.london, exprs.unique)
degs <- filterDEGS(degs, 0.05, 0.65, adj=F)
