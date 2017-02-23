source("degs_utils.R")
source("plots_utils.R")
library(cowplot)
library(illuminaHumanv3.db)
library(psych)

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
exprs <- read.table("../exprs/exprs_all.tsv", sep="\t", quote = "",
                    header=TRUE, check.names=FALSE)
pdata  <- read.table("../pdata/pdata.tsv", sep="\t", quote = "",
                   header=TRUE, check.names=FALSE)

## Filter lists by p-val
names <- LC.degs[,c(1:3)]
LC.degs <- filterDEGS(LC.degs, 0.01, 0)[,c(1, 4)]
LP.degs <- filterDEGS(LP.degs, 0.01, 0)[,c(1, 4)]
HC.degs <- filterDEGS(HC.degs, 0.01, 0)[,c(1, 4)]
HP.degs <- filterDEGS(HP.degs, 0.01, 0)[,c(1, 4)]
CP.degs <- filterDEGS(CP.degs, 0.05, 0)[,c(1, 4)]
colnames(LC.degs) <- c("ENTREZID", "logFC.LC")
colnames(LP.degs) <- c("ENTREZID", "logFC.LP")
colnames(HC.degs) <- c("ENTREZID", "logFC.HC")
colnames(HP.degs) <- c("ENTREZID", "logFC.HP")
colnames(CP.degs) <- c("ENTREZID", "logFC.CP")

## Merge all lists
merged <- merge(names, LC.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, LP.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, HC.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, HP.degs, by="ENTREZID", all=TRUE)
merged <- merge(merged, CP.degs, by="ENTREZID", all=TRUE)

## Filter if all logFC is NA
ind <- which(!rowSums(!is.na(merged[,4:ncol(merged)]))) 
merged <- merged[-ind,]

## Add average exprs per group per gene
exprs <- exprs[which(rownames(exprs) %in% merged$ENTREZID),]
exprs <- exprs[order(match(rownames(exprs), merged$ENTREZID)),]

p <- pdata[pdata$Condition=="Low risk",]
e <- rowMeans(exprs[,which(colnames(exprs) %in% rownames(p))])
merged$mL <- e

p <- pdata[pdata$Condition=="High risk",]
e <- rowMeans(exprs[,which(colnames(exprs) %in% rownames(p))])
merged$mH <- e

p <- pdata[pdata$Condition=="Control",]
e <- rowMeans(exprs[,which(colnames(exprs) %in% rownames(p))])
merged$mC <- e

p <- pdata[pdata$Condition=="Preeclampsia",]
e <- rowMeans(exprs[,which(colnames(exprs) %in% rownames(p))])
merged$mP <- e

## Find the difference between logFCs
merged$logFC.LC_LP <- round(merged$logFC.LC - merged$logFC.LP, 4)
merged$logFC.LC_HC <- round(merged$logFC.LC - merged$logFC.HC, 4)
merged$logFC.LC_HP <- round(merged$logFC.LC - merged$logFC.HP, 4)
merged$logFC.LP_HP <- round(merged$logFC.LP - merged$logFC.HP, 4)
merged$logFC.HC_HP <- round(merged$logFC.HC - merged$logFC.HP, 4)

mm <- merged
merged <- mm
## Filter alt
total <- filterJointDEGS(merged)
tfs <- total[which(total$ENTREZID %in% TFs$ENTREZID),]
total$isTF <- total$ENTREZID %in% TFs$ENTREZID
total <- total[,c(1:3, 18, 4:17)]
write.table(total, "../degs/interesting_genes.tsv", sep="\t", row.names = FALSE, quote=FALSE)

############### GRN
total <- filterJointDEGS(merged, 0.7, 0.8, 0.8)
total <- total[,c(1:4)]
tfs <- as.character(total[which(total$ENTREZID %in% TFs$ENTREZID),]$SYMBOL)
LC.degs <- read.table("../degs/LC_degs.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)
LC.degs <- LC.degs[which(LC.degs$ENTREZID %in% total$ENTREZID), ]
LC.degs <- LC.degs[, -c(1,3:9)]

HP.degs <- read.table("../degs/HP_degs.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)
HP.degs <- HP.degs[which(HP.degs$ENTREZID %in% total$ENTREZID), ]
HP.degs <- HP.degs[, -c(1,3:9)]

write.table(t(tfs), "../grn/regulators.txt", sep=" ", col.names=FALSE, row.names = FALSE, quote=FALSE)
write.table(LC.degs, "../grn/norma.tsv", sep="\t", row.names = FALSE, quote=FALSE)
write.table(HP.degs, "../grn/disease.tsv", sep="\t", row.names = FALSE, quote=FALSE)
###############

### Check LH tfs

LH.degs <- read.table("../degs/LH_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)
LH.degs <- filterDEGS(LH.degs, 0.05, 0.7, adj=FALSE)
tfs <- as.character(LH.degs[which(LH.degs$ENTREZID %in% TFs$ENTREZID),]$SYMBOL)

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

### Barplot for particular gene
exprs <- read.table("../exprs/exprs_all.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)
pdata  <- read.table("../pdata/pdata.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)

symbol = "FLT1"
pl <- geneBarPlot(exprs, pdata, symbol)
save_plot(paste("../plots/GenesBP/", symbol, ".pdf", sep=""),
          base_height=5, base_aspect_ratio=1.5, pl)
