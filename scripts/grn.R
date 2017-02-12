library(GGally)
library(network)
library(sna)
library(ggplot2)
library( RColorBrewer)

### Read bnf output
disease <- read.table("../grn/disease.sif", sep="\t", quote = "", header=FALSE, check.names=FALSE)
norma <- read.table("../grn/norma.sif", sep="\t", quote = "", header=FALSE, check.names=FALSE)
tfs <- as.character(read.table("../grn/regulators.txt", sep=" ", quote = "", header=FALSE, check.names=FALSE,
                  stringsAsFactors = FALSE))

colnames(disease) <- c("IN", "weight", "OUT")
colnames(norma) <- c("IN", "weight", "OUT")

CP.degs <- read.table("../degs/CP_degs.tsv", sep="\t", quote = "",
                      header=TRUE, check.names=FALSE)


### Filter

CP.degs <- filterDEGS(CP.degs, 0.05, 0.5)

cutoff = 0.7
disease <- disease[order(disease$weight, decreasing = TRUE),]
disease <- disease[which(disease$weight>cutoff),]
disease <- disease[, c(1, 3, 2)]

norma <- norma[order(norma$weight, decreasing = TRUE),]
norma <- norma[which(norma$weight>cutoff),]
norma <- norma[, c(1, 3, 2)]

# Intersect edges
norma$common <- c('FALSE', 'TRUE')[(paste(norma$IN, norma$OUT) %in% paste(disease$IN, disease$OUT))+1]
disease$common <- c('FALSE', 'TRUE')[paste(disease$IN, disease$OUT) %in% (paste(norma$IN, norma$OUT))+1]


# network
d <- network(disease, directed=TRUE, matrix.type="edgelist", ignore.eval = FALSE)
isTF <- get.vertex.attribute(d,"vertex.names") %in% tfs
isCP <- get.vertex.attribute(d,"vertex.names") %in% CP.degs$SYMBOL
color <- c()
for (i in 1:length(isTF)) {
  if ((isTF[i]==FALSE)&&(isCP[i]==FALSE)) {color <- c(color, "steelblue")}
  if ((isTF[i]==TRUE)&&(isCP[i]==FALSE)) {color <- c(color, "tomato")}
  if ((isTF[i]==FALSE)&&(isCP[i]==TRUE)) {color <- c(color, "mediumseagreen")}
  if ((isTF[i]==TRUE)&&(isCP[i]==TRUE)) {color <- c(color, "#DC143C")}
}
d %v% "group" = color
set.edge.attribute(d, "color", ifelse(d %e% "common" == FALSE, "gray", "tomato"))

n <- network(norma, directed=TRUE, matrix.type="edgelist", ignore.eval = FALSE)
isTF <- get.vertex.attribute(n,"vertex.names") %in% tfs
isCP <- get.vertex.attribute(n,"vertex.names") %in% CP.degs$SYMBOL
color <- c()
for (i in 1:length(isTF)) {
  if ((isTF[i]==FALSE)&&(isCP[i]==FALSE)) {color <- c(color, "steelblue")}
  if ((isTF[i]==TRUE)&&(isCP[i]==FALSE)) {color <- c(color, "tomato")}
  if ((isTF[i]==FALSE)&&(isCP[i]==TRUE)) {color <- c(color, "mediumseagreen")}
  if ((isTF[i]==TRUE)&&(isCP[i]==TRUE)) {color <- c(color, "#DC143C")}
}
n %v% "color" = color
set.edge.attribute(n, "color", ifelse(n %e% "common" == FALSE, "gray", "tomato"))


pl <- ggnet2(d, node.size = 16, color = "gray15", edge.color = "color", node.alpha = 0.7,
             label = TRUE, label.size = 3, label.color="group") +
      theme(panel.background = element_rect(fill = "grey15"))
save_plot("../plots/GRN/disease.pdf",
          base_height=7, base_aspect_ratio=1.3, pl)

pl <- ggnet2(n, node.size = 16, color = "gray15", edge.color = "color", node.alpha = 0.7,
             label = TRUE, label.size = 3, label.color="color") +
      theme(panel.background = element_rect(fill = "grey15"))
save_plot("../plots/GRN/norma.pdf",
          base_height=7, base_aspect_ratio=1.3, pl)
