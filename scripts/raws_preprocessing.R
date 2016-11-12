library(cowplot)
library(ggfortify)
library(beadarray)
library(arrayQualityMetrics)
library(sva)
source("pcaPlots.R")

# General
studies <- read.table("../pdata/studies.tsv", header = TRUE, sep = "\t")

# London
i=1
path = paste("../raws/", studies[i,]$ID, sep="")
pdata = read.table(paste("../pdata/", studies[i,]$ID, "_pdata.tsv", sep=""), sep="\t", head=TRUE, stringsAsFactors = FALSE)
idatfiles = dir(path, pattern="idat", full.name=TRUE)
pdata <- pdata[c(-23, -24),]
idatfiles <- idatfiles[c(-23, -24)]

susp = which(pdata$Comment=="Suspicious")
pdata <- pdata[-susp,]
idatfiles <- idatfiles[-susp]

raw.data <- readIdatFiles(idatfiles)
raw.data@phenoData = AnnotatedDataFrame(pdata)
N.data <- normaliseIllumina(channel(raw.data, "Green"), method="neqc", transform="none")

# Batch-effect removal
batch = as.factor(pData(N.data)$Sentrix_ID)
mod = model.matrix(~as.factor(Condition), data=pData(N.data))
#mod = model.matrix(~1, data=pData(N.data))
combat_edata = ComBat(dat=exprs(N.data), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
N.data@assayData$exprs <- combat_edata

# Perform PCA and plot
pca = prcomp(t(exprs(N.data)))
pl <- pcaPlots(pca, pData(N.data), c("Condition", "Sentrix_ID"), studies[i,]$ID)

# Save plot for manual quality control
save_plot(paste("../plots/", studies[i,]$ID, "_PCA.pdf", sep=""),
          base_height=4, base_aspect_ratio = 2.2, pl)

# AQM
arrayQualityMetrics(expressionset = N.data,
                    outdir = "../plots/Report_1",
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Condition"))
