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
LC.degs <- read.table("../degs/LC_degs.tsv", sep="\t", quote = "",
                             header=TRUE, check.names=FALSE)
HP.degs <- read.table("../degs/HP_degs.tsv", sep="\t", quote = "",
                             header=TRUE, check.names=FALSE)

## Filter and order lists
LC.degs <- filterDEGS(LC.degs, 0.01, 0)
HP.degs <- filterDEGS(HP.degs, 0.01, 0)

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
