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
LC.degs <- read.table("../degs/LC_degs.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)
HP.degs <- read.table("../degs/HP_degs.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)
CP.degs <- read.table("../degs/CP_degs.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)
LH.degs <- read.table("../degs/LH_degs.tsv", sep="\t", quote = "", header=TRUE, check.names=FALSE)

## Filter and order lists
LC.degs <- filterDEGS(LC.degs, 0.01, 0)
HP.degs <- filterDEGS(HP.degs, 0.01, 0)
#CP.degs <- filterDEGS(CP.degs, 0.05, 0)
#LH.degs <- filterDEGS(LH.degs, 0.05, 0, adj = FALSE)

### GSEA GO analysis

## LC and HP
lc.degs.gse <- LC.degs$logFC
names(lc.degs.gse) <- LC.degs$ENTREZID
lc.degs.gse <- lc.degs.gse[order(lc.degs.gse, decreasing = TRUE)]

hp.degs.gse <- HP.degs$logFC
names(hp.degs.gse) <- HP.degs$ENTREZID
hp.degs.gse <- hp.degs.gse[order(hp.degs.gse, decreasing = TRUE)]

doGSEGO(lc.degs.gse, hp.degs.gse, by="GeneRatio", order=TRUE)

### EnrichGO

## CP
cp.degs.ego <- filterDEGS(CP.degs, 0.05, 0.5)
cp.degs.ego <- cp.degs.ego$ENTREZID

cp.ego <- enrichGO(cp.degs.ego, universe = as.character(CP.degs$ENTREZID), OrgDb = org.Hs.eg.db, ont = "CC")
cp.ego.filtered <- dropGO(cp.ego, level=1)
cp.ego.filtered <- dropGO(cp.ego.filtered, level=2)
pl <- dotplot(cp.ego.filtered, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/CP_eGO_CC.pdf", base_height=5, base_aspect_ratio=2.1, pl)
save_plot("../plots/FunctionalAnalysis/CP_eGO_CC.svg", base_height=5, base_aspect_ratio=2.1, pl)

cp.ego <- enrichGO(cp.degs.ego, universe = as.character(CP.degs$ENTREZID), OrgDb = org.Hs.eg.db, ont = "BP")
cp.ego.filtered <- dropGO(cp.ego, level=1)
cp.ego.filtered <- dropGO(cp.ego.filtered, level=2)
pl <- dotplot(cp.ego.filtered, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/CP_eGO_BP.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/CP_eGO_BP.svg", base_height=5, base_aspect_ratio=1.7, pl)

cp.ego <- enrichGO(cp.degs.ego, universe = as.character(CP.degs$ENTREZID), OrgDb = org.Hs.eg.db, ont = "MF")
cp.ego.filtered <- dropGO(cp.ego, level=1)
cp.ego.filtered <- dropGO(cp.ego.filtered, level=2)
pl <- dotplot(cp.ego.filtered, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/CP_eGO_MF.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/CP_eGO_MF.svg", base_height=5, base_aspect_ratio=1.7, pl)

## LH
lh.degs.ego <- filterDEGS(LH.degs, 0.05, 0.5, adj = FALSE)
lh.degs.ego <- lh.degs.ego$ENTREZID

lh.ego <- enrichGO(lh.degs.ego, universe = as.character(LH.degs$ENTREZID), OrgDb = org.Hs.eg.db, ont = "CC")
lh.ego.filtered <- dropGO(lh.ego, level=1)
lh.ego.filtered <- dropGO(lh.ego.filtered, level=2)
pl <- dotplot(lh.ego, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/LH_eGO_CC.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/LH_eGO_CC.svg", base_height=5, base_aspect_ratio=1.7, pl)

lh.ego <- enrichGO(lh.degs.ego, universe = as.character(LH.degs$ENTREZID), OrgDb = org.Hs.eg.db, ont = "BP")
lh.ego.filtered <- dropGO(lh.ego, level=1)
lh.ego.filtered <- dropGO(lh.ego.filtered, level=2)
pl <- dotplot(lh.ego.filtered, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/LH_eGO_BP.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/LH_eGO_BP.svg", base_height=5, base_aspect_ratio=1.7, pl)

lh.ego <- enrichGO(lh.degs.ego, universe = as.character(LH.degs$ENTREZID), OrgDb = org.Hs.eg.db, ont = "MF")
lh.ego.filtered <- dropGO(lh.ego, level=1)
lh.ego.filtered <- dropGO(lh.ego.filtered, level=2)
pl <- dotplot(lh.ego.filtered, showCategory = 15, font.size=10)
save_plot("../plots/FunctionalAnalysis/LH_eGO_MF.pdf", base_height=5, base_aspect_ratio=2.1, pl)
save_plot("../plots/FunctionalAnalysis/LH_eGO_MF.svg", base_height=5, base_aspect_ratio=2.1, pl)

### EnrichKEGG

## CP
cp.kk <- enrichKEGG(gene = cp.degs.ego)
pl <- dotplot(cp.kk, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/CP_kk.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/CP_kk.svg", base_height=5, base_aspect_ratio=1.7, pl)

## LH
lh.kk <- enrichKEGG(gene = lh.degs.ego)
pl <- dotplot(lh.kk, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/LH_kk.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/LH_kk.svg", base_height=5, base_aspect_ratio=1.7, pl)

### EnrichPathway

## CP
cp.path <- enrichPathway(gene = cp.degs.ego)
pl <- dotplot(cp.path, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/CP_path.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/CP_path.svg", base_height=5, base_aspect_ratio=1.7, pl)

## LH
lh.path <- enrichPathway(gene = lh.degs.ego)
pl <- dotplot(lh.path, showCategory = 15)
save_plot("../plots/FunctionalAnalysis/LH_path.pdf", base_height=5, base_aspect_ratio=1.7, pl)
save_plot("../plots/FunctionalAnalysis/LH_path.svg", base_height=5, base_aspect_ratio=1.7, pl)

### CompareClusters

LC.degs.filtered <- filterDEGS(LC.degs, 0.01, 0.7)
HP.degs.filtered <- filterDEGS(HP.degs, 0.01, 0.7)
gc <- list(LC=as.character(LC.degs.filtered$ENTREZID), HP=as.character(HP.degs.filtered$ENTREZID))

## enrichPathway
ck <- compareCluster(geneCluster = gc, fun = "enrichPathway")
pl <- dotplot(ck, showCategory = 30) + theme_bw(base_size = 18) +
      theme(axis.title.x=element_text(margin=margin(t=10)))

save_plot("../plots/FunctionalAnalysis/LCHP_enrichPathway.pdf",
          base_height=20, base_aspect_ratio=0.7, pl)
save_plot("../plots/FunctionalAnalysis/LCHP_enrichPathway.svg",
          base_height=20, base_aspect_ratio=0.7, pl)

## enrichGO
ck <- compareCluster(geneCluster = gc, fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "CC")
pl <- dotplot(ck, showCategory = 50) + theme_bw(base_size = 18) +
  theme(axis.title.x=element_text(margin=margin(t=10)))

save_plot("../plots/FunctionalAnalysis/LCHP_enrichGO_CC.pdf",
          base_height=15, base_aspect_ratio=0.7, pl)
save_plot("../plots/FunctionalAnalysis/LCHP_enrichGO_CC.svg",
          base_height=15, base_aspect_ratio=0.7, pl)
