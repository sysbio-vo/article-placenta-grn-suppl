library(org.Hs.eg.db)
library(stringr)
library(cowplot)
library(reshape2)
library(RColorBrewer)
source("degs_utils.R")
source("go_utils.R")
library(devtools)
library(ReactomePA)
#dev_mode(on=T)
#install_github("GuangchuangYu/DOSE")
#dev_mode(on=F)
library(DOSE)
library(clusterProfiler)

### Reading DEG lists
LC.degs.unique <- read.table("../degs/LC_degs_absolute.tsv", sep="\t",  quote = "",
                             header=TRUE, check.names=FALSE)
HP.degs.unique <- read.table("../degs/HP_degs_absolute.tsv", sep="\t",  quote = "",
                             header=TRUE, check.names=FALSE)
LC.degs <- read.table("../degs/LC_degs.tsv", sep="\t", quote = "",
                             header=TRUE, check.names=FALSE)
HP.degs <- read.table("../degs/HP_degs.tsv", sep="\t", quote = "",
                             header=TRUE, check.names=FALSE)
CP.degs <- read.table("../degs/CP_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)

## Filter and order lists
LC.degs <- filterDEGS(LC.degs, 0.01, 0)
LC.degs <- LC.degs[order(LC.degs$logFC, decreasing = TRUE),]

HP.degs <- filterDEGS(HP.degs, 0.01, 0)
HP.degs <- HP.degs[order(HP.degs$logFC, decreasing = TRUE),]

LC.degs.unique <- LC.degs.unique[order(LC.degs.unique$logFC, decreasing = TRUE),]
HP.degs.unique <- HP.degs.unique[order(HP.degs.unique$logFC, decreasing = TRUE),]

### GSEA GO analysis

lc.degs.gse <- LC.degs$logFC
names(lc.degs.gse) <- LC.degs$ENTREZID
hp.degs.gse <- HP.degs$logFC
names(hp.degs.gse) <- HP.degs$ENTREZID

doGSEGO(lc.degs.gse, hp.degs.gse, by="GeneRatio", order=TRUE)

### CompareClusters

LC.degs.filtered <- filterDEGS(LC.degs, 0.01, 0.7)
HP.degs.filtered <- filterDEGS(HP.degs, 0.01, 0.7)
gc <- list(LC=as.character(LC.degs.filtered$ENTREZID), HP=as.character(HP.degs.filtered$ENTREZID))

## enrichPathway
ck <- compareCluster(geneCluster = gc, fun = "enrichPathway")
pl <- dotplot(ck, showCategory = 30) + theme_bw(base_size = 18) +
      theme(axis.title.x=element_text(margin=margin(t=10)))

save_plot("../plots/FunctionalAnalysis/LCPH_enrichPathway.pdf",
          base_height=20, base_aspect_ratio=0.7, pl)
save_plot("../plots/FunctionalAnalysis/LCPH_enrichPathway.svg",
          base_height=20, base_aspect_ratio=0.7, pl)

## enrichGO
ck <- compareCluster(geneCluster = gc, fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "CC")
pl <- dotplot(ck, showCategory = 50) + theme_bw(base_size = 18) +
  theme(axis.title.x=element_text(margin=margin(t=10)))

save_plot("../plots/FunctionalAnalysis/LCPH_enrichGO_CC.pdf",
          base_height=15, base_aspect_ratio=0.7, pl)
save_plot("../plots/FunctionalAnalysis/LCPH_enrichGO_CC.svg",
          base_height=15, base_aspect_ratio=0.7, pl)

### Compare LC and HP

# Find common genes for LC and HP
common <- intersect(LC.degs$SYMBOL, HP.degs$SYMBOL)
LC.degs.c <- LC.degs[(LC.degs$SYMBOL %in% common), c("SYMBOL", "GENENAME", "logFC", "AveExpr")]
HP.degs.c <- HP.degs[(HP.degs$SYMBOL %in% common), c("SYMBOL", "GENENAME", "logFC", "AveExpr")]

# Merge lists
merged <- merge(LC.degs.c, HP.degs.c, by="SYMBOL")
merged <- merged[, -5]
colnames(merged) <- c("SYMBOL", "GENENAME", "logFC.LC", "AveExpr.LC", "logFC.HP", "AveExpr.HP")

# Check correlation
cor(merged$logFC.LC, merged$logFC.HP)

# Find difference between LC and HP logFC and filter based on difference
merged$logFC.HP_LC <- -(merged$logFC.LC-merged$logFC.HP)
merged.filtered <- merged[abs(merged$logFC.HP_LC)>1,]

# Search for the same genes in CP DEG list
CP.degs.f <- CP.degs[CP.degs$SYMBOL %in% merged.filtered$SYMBOL,]
CP.degs.f <- CP.degs.f[order(match(CP.degs.f$SYMBOL, merged.filtered$SYMBOL)),c("logFC", "adj.P.Val")]
colnames(CP.degs.f)[1] <- c("logFC.CP")
CP.degs.f$logFC.CP[which(CP.degs.f$adj.P.Val >= 0.05)] <- NA

merged.filtered <- cbind(merged.filtered, CP.degs.f)

## Prepare for plotting
m <- merged.filtered[,c("SYMBOL", "logFC.LC", "logFC.HP", "logFC.HP_LC", "logFC.CP")]
rownames(m) <- merged.filtered$SYMBOL
m <- melt(m, id.vars="SYMBOL")
m$value <- round(m$value, 2)
# Create separate variable to use with facetting
m$LCHP = "LC minus HP logFC compared to CP"
ind <- c(which(m$variable=="logFC.LC"), which(m$variable=="logFC.HP"))
m$LCHP[ind] <- "LC and HP logFC"
# Sort by logFC.LC
sort <- m[m$variable=="logFC.LC",]
sort <- sort[order(sort$value, decreasing = TRUE),]
m$SYMBOL <- factor(m$SYMBOL, levels=sort$SYMBOL)

# Use custom palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

pl <- ggplot(m, aes(x=SYMBOL, y=variable, fill=value)) + coord_flip() +
      geom_tile(aes(fill = value)) + 
      geom_text(aes(label = value)) +
      scale_fill_gradientn(colours = myPalette(100), na.value="white") +
      theme(axis.title = element_blank()) +
      labs(fill='logFC') +
      facet_grid(~LCHP, scales = "free_x") +
      theme(axis.text.y = element_text(size=10),
            panel.spacing.x = unit(5, "mm")) 

save_plot("../plots/FunctionalAnalysis/LCPH_heatmap.pdf",
          base_height=14, base_aspect_ratio=0.57, pl)
save_plot("../plots/FunctionalAnalysis/LCPH_heatmap.svg",
          base_height=14, base_aspect_ratio=0.57, pl)
