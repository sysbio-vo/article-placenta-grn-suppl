str_count <- function(string, pattern="") {
  sapply(string, str_count_item, pattern=pattern)
}

str_count_item <- function(string, pattern = "") {
  length(gregexpr(pattern, string)[[1]])
}

fortify.gsea <- function(data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, split=NULL) {
  res <- as.data.frame(data)
  res$Count <- str_count(res$core_enrichment, "/") + 1
  res$.sign <- "activated"
  res$.sign[res$NES < 0] <- "suppressed"
  if (drop) {
    res <- res[res$Count != 0, ]
  }
  
  res$GeneRatio <- res$Count / res$setSize
  
  if (order) {
    if (by == "Count") {
      idx <- order(res$Count, decreasing=TRUE)
    } else {
      idx <- order(res$GeneRatio, decreasing=TRUE)
    }
    res <- res[idx,]
  }
  
  topN <- function(res, showCategory) {
    if ( is.numeric(showCategory) ) {
      if ( showCategory <= nrow(res) ) {
        res <- res[1:showCategory,]
      }
    } else { ## selected categories
      res <- res[res$ID %in% showCategory,]
    }
    return(res)
  }
  
  if (is.null(split)) {
    res <- topN(res, showCategory)
  } else {
    lres <- split(res, as.character(res[, split]))
    lres <- lapply(lres, topN, showCategory = showCategory)
    res <- do.call('rbind', lres)
  }
  
  res$Description <- factor(res$Description,
                            levels=rev(res$Description))
  
  return(res)
}

dotplot.gsea <- function(df, font.size=12) {
  require(stringr)
  idx <- order(df$GeneRatio, decreasing = FALSE)
  df$Description <- str_wrap(df$Description, width = 70)
  d <- df$Description[idx]
  d <- d[!duplicated(d)]
  df$Description <- factor(df$Description, levels = d)

  x <- "GeneRatio"
  size <- "Count"
  colorBy="p.adjust"

  pl <- ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() + scale_color_gradient(low="red", high="blue") +
    ylab("") + facet_grid(.~.sign + Group, scales="fixed") + theme_bw(base_size = 18) +
    theme(axis.title.x=element_text(margin=margin(t=10)),
          axis.text.y = element_text(lineheight=0.75),
          panel.spacing.x = unit(5, "mm"))  
  return(pl)
}

doGSEGO <- function(x, y, x.group="LC", y.group="HP", showCategory=20, by="Count", order=TRUE) {
  go <- c("BP", "CC", "MF")
  for (i in go) {
    x.gse <- gseGO(x, OrgDb = org.Hs.eg.db, ont = i)
    y.gse <- gseGO(y, OrgDb = org.Hs.eg.db, ont = i)
    if (nrow(x.gse)>0 && nrow(y.gse)>0) {
      x.df <- fortify.gsea(x.gse, showCategory = 20, order=order, split=".sign", by=by)
      y.df <- fortify.gsea(y.gse, showCategory = 20, order=order, split=".sign", by=by)
      x.df$Group <- rep(x.group, nrow(x.df))
      y.df$Group <- rep(y.group, nrow(y.df))
      df <- rbind(x.df, y.df)
      pl <- dotplot.gsea(df)
      save_plot(paste("../plots/FunctionalAnalysis/", x.group, y.group, "_gseGO_", i, "_", by, ".pdf", sep=""),
                base_height=25, base_aspect_ratio=0.7, pl)
      save_plot(paste("../plots/FunctionalAnalysis/", x.group, y.group, "_gseGO_", i, "_", by, ".svg", sep=""),
                base_height=25, base_aspect_ratio=0.7, pl)
    }
  }
}