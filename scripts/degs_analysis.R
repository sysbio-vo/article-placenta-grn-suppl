source("degs_utils.R")
source("plots_utils.R")
library(cowplot)


### Reading DEG lists
LC.degs <- read.table("../degs/LC_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)
LP.degs <- read.table("../degs/LP_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)
HC.degs <- read.table("../degs/HC_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)
HP.degs <- read.table("../degs/HP_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)
CP.degs <- read.table("../degs/CP_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)
TFs  <- read.table("../pdata/TFs.tsv", sep="\t", quote = "",
                   header=TRUE, check.names=FALSE)

## Filter lists by p-val
names <- LC.degs[,c(1:3)]
LC.degs <- filterDEGS(LC.degs, 0.01, 0)[,c(1, 4:5)]
LP.degs <- filterDEGS(LP.degs, 0.01, 0)[,c(1, 4:5)]
HC.degs <- filterDEGS(HC.degs, 0.01, 0)[,c(1, 4:5)]
HP.degs <- filterDEGS(HP.degs, 0.01, 0)[,c(1, 4:5)]
CP.degs <- filterDEGS(CP.degs, 0.05, 0)[,c(1, 4:5)]
colnames(LC.degs) <- c("ENTREZID", "logFC.LC", "AveExprs.LC")
colnames(LP.degs) <- c("ENTREZID", "logFC.LP", "AveExprs.LP")
colnames(HC.degs) <- c("ENTREZID", "logFC.HC", "AveExprs.HC")
colnames(HP.degs) <- c("ENTREZID", "logFC.HP", "AveExprs.HP")
colnames(CP.degs) <- c("ENTREZID", "logFC.CP", "AveExprs.CP")

## Merge all lists
merged <- merge(names, LC.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, LP.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, HC.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, HP.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, CP.degs, by="ENTREZID", all=TRUE)

## Filter if all logFC is NA
ind <- which(!rowSums(!is.na(merged[,4:ncol(merged)]))) 
merged <- merged[-ind,]
merged <- merged[,c(1:3, 4, 6, 8, 10, 12, 5, 7, 9, 11, 13)]

## Find the difference between logFCs
merged$logFC.LC_LP <- merged$logFC.LC - merged$logFC.LP
merged$logFC.LC_HC <- merged$logFC.LC - merged$logFC.HC
merged$logFC.LC_HP <- merged$logFC.LC - merged$logFC.HP
merged$logFC.LP_HC <- merged$logFC.LP - merged$logFC.HC
merged$logFC.LP_HP <- merged$logFC.LP - merged$logFC.HP
merged$logFC.HC_HP <- merged$logFC.HC - merged$logFC.HP

## Filter alt
tfs.big <- merged[which(merged$ENTREZID %in% TFs$ENTREZID),]
ind <- which(apply(merged[, c(14:19)], MARGIN = 1, function(x) any(abs(x) > 1, na.rm = TRUE))==TRUE)
merged.diff <- merged[ind,]
merged <- merged[-ind,]

ind <- which(apply(merged[, c(4:7)], MARGIN = 1, function(x) all(abs(x) > 0, na.rm = FALSE))==TRUE)
merged.na <- merged[-ind,]
merged <- merged[ind,]

ind <- which(apply(merged.na[, c(4:8)], MARGIN = 1, function(x) all(abs(x) < 0.7, na.rm = TRUE))==TRUE)
merged.na <- merged.na[-ind,]

merged <- merged[which(abs(merged$logFC.CP)>0.7),]

total <- rbind(merged.diff, merged.na, merged)
tfs <- total[which(total$ENTREZID %in% TFs$ENTREZID),]

### Compare LC and HP

# Find common genes for LC and HP
LCHP <- merged[,c("SYMBOL", "logFC.LC", "logFC.HP", "logFC.LC_HP", "logFC.CP")]
LCHP <- LCHP[!is.na(rowSums(LCHP[,2:3])),]

# Check correlation
cor(LCHP$logFC.LC, LCHP$logFC.HP)

# Filter based on difference between LC and HP
LCHP <- LCHP[abs(LCHP$logFC.LC_HP)>1,]
LCHP$logFC.LC_HP <- -LCHP$logFC.LC_HP
colnames(LCHP)[4] <- "logFC.HP_LC"

## Plot heatmap
pl <- genesHeatmap(LCHP, "LC", "HP")

save_plot("../plots/FunctionalAnalysis/LCPH_heatmap.pdf",
          base_height=14, base_aspect_ratio=0.57, pl)
save_plot("../plots/FunctionalAnalysis/LCPH_heatmap.svg",
          base_height=14, base_aspect_ratio=0.57, pl)


### Compare LC and LP

# Find common genes for LC and LP
LCLP <- merged[,c("SYMBOL", "logFC.LC", "logFC.LP", "logFC.LC_LP", "logFC.CP")]
LCLP <- LCLP[!is.na(rowSums(LCLP[,2:3])),]

# Check correlation
cor(LCLP$logFC.LC, LCLP$logFC.LP)

# Filter based on difference between LC and HP
LCLP <- LCLP[abs(LCLP$logFC.LC_LP)>0.82,]
LCLP$logFC.LC_LP <- -LCLP$logFC.LC_LP
colnames(LCLP)[4] <- "logFC.LP_LC"

## Plot heatmap
pl <- genesHeatmap(LCLP, "LC", "LP")

save_plot("../plots/FunctionalAnalysis/LCLP_heatmap.pdf",
          base_height=14, base_aspect_ratio=0.57, pl)
save_plot("../plots/FunctionalAnalysis/LCLP_heatmap.svg",
          base_height=14, base_aspect_ratio=0.57, pl)
