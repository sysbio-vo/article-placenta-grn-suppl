require(ggplot2)
require(cowplot)
require(ggfortify)

pcaPlots <- function(pca.data, pheno.data, meta.vars, title) {
  pheno.data[] <- lapply(pheno.data, as.character)
  title <- ggdraw() + draw_label(title, fontface='bold')  
  plots <- c()
  for (i in meta.vars) {
    pl <- autoplot(pca.data, data = pheno.data, colour=i)
    plots <- c(plots, list(pl))
  }
  pl <- plot_grid(plotlist = plots, ncol=length(meta.vars))
  pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  return(pl)
}