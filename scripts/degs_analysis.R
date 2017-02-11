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
merged$meanL <- e

p <- pdata[pdata$Condition=="High risk",]
e <- rowMeans(exprs[,which(colnames(exprs) %in% rownames(p))])
merged$meanH <- e

p <- pdata[pdata$Condition=="Control",]
e <- rowMeans(exprs[,which(colnames(exprs) %in% rownames(p))])
merged$meanC <- e

p <- pdata[pdata$Condition=="Preeclampsia",]
e <- rowMeans(exprs[,which(colnames(exprs) %in% rownames(p))])
merged$meanP <- e

## Find the difference between logFCs
merged$logFC.LC_LP <- round(merged$logFC.LC - merged$logFC.LP, 4)
merged$logFC.LC_HC <- round(merged$logFC.LC - merged$logFC.HC, 4)
merged$logFC.LC_HP <- round(merged$logFC.LC - merged$logFC.HP, 4)
merged$logFC.LP_HP <- round(merged$logFC.LP - merged$logFC.HP, 4)
merged$logFC.HC_HP <- round(merged$logFC.HC - merged$logFC.HP, 4)

mm <- merged
merged <- mm
## Filter alt
tfs.big <- merged[which(merged$ENTREZID %in% TFs$ENTREZID),]
ind <- which(apply(merged[, c(13:17)], MARGIN = 1, function(x) any(abs(x) > 1, na.rm = TRUE))==TRUE)
merged.diff <- merged[ind,]
merged <- merged[-ind,]

ind <- which(apply(merged[, c(4:7)], MARGIN = 1, function(x) all(abs(x) > 0, na.rm = FALSE))==TRUE)
merged.na <- merged[-ind,]
merged <- merged[ind,]

ind <- which(apply(merged.na[, c(4:8)], MARGIN = 1, function(x) all(abs(x) < 1, na.rm = TRUE))==TRUE)
merged.na <- merged.na[-ind,]

merged <- merged[which(abs(merged$logFC.CP)>1),]

total <- rbind(merged.diff, merged.na, merged)
tfs <- total[which(total$ENTREZID %in% TFs$ENTREZID),]
total$isTF <- total$ENTREZID %in% TFs$ENTREZID
total <- total[,c(1:3, 18, 4:17)]

write.table(total, "../degs/interesting_genes.tsv", sep="\t", row.names = FALSE, quote=FALSE)

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
