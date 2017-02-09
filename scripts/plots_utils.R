require(ggplot2)
require(cowplot)
require(ggfortify)
library(grid)
library(reshape2)
library(RColorBrewer)

getAspectRatio <- function(p){
  gb <- ggplot_build(p)
  g <- ggplot_gtable(gb)
  
  nullw <- sapply(g$widths, attr, "unit")
  nullh <- sapply(g$heights, attr, "unit")
  
  # ar of plot
  if(any(nullw == "null"))
    ar <- unlist(g$widths[nullw == "null"]) / unlist(g$heights[nullh == "null"])

  # ar of plot + legend
  g$fullwidth <- convertWidth(sum(g$widths), "in", valueOnly=TRUE)
  g$fullheight <- convertHeight(sum(g$heights), "in", valueOnly=TRUE)
  ar <- g$fullwidth / g$fullheight

  return(ar)
}

pcaPlots <- function(pca.data, pheno.data, meta.vars, title, ncol) {
  pheno.data[] <- lapply(pheno.data, as.character)
  plots <- c()
  ar <- -100
  for (i in meta.vars) {
    pl <- autoplot(pca.data, data = pheno.data, colour=i) +
                  coord_fixed()
    newar <- getAspectRatio(pl)
    if (ar<newar) {
      ar <- newar
    }
    plots <- c(plots, list(pl))
  }
  if(missing(ncol)) {
    pl <- plot_grid(plotlist = plots, ncol=length(meta.vars), align="hv")
  } else {
    pl <- plot_grid(plotlist = plots, ncol=ncol, align="hv")
  }
  if (!missing(title)) {
    title <- ggdraw() + draw_label(title, fontface='bold')  
    pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  }
  pl <- pl + theme(plot.margin=margin(t=10, r=10, b=10, l=10))
  return(list(pl, ar))
}

genesHeatmap <- function(df, nameA, nameB) {
  rownames(df) <- df$SYMBOL
  m <- melt(df, id.vars="SYMBOL")
  m$value <- round(m$value, 2)
  # Create separate variable to use with facetting
  m$LCHP = paste("2. ", nameB, " minus ", nameA, " logFC compared to CP", sep="")
  ind <- c(which(m$variable==paste("logFC.", nameA, sep="")),
           which(m$variable==paste("logFC.", nameB, sep="")))
  m$LCHP[ind] <- paste("1. ", nameA, " and ", nameB, " logFC", sep="")
  # Sort by logFC.LC
  sort <- m[m$variable==paste("logFC.", nameA, sep=""),]
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
          panel.spacing.x = unit(7, "mm"))   
  
  return(pl)
}