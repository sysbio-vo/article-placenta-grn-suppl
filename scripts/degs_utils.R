getDEGS <- function(meta.vars, pheno.data, exprs, noannotate=FALSE) {
  # Subset groups for comparison
  ind <- c()
  for (i in meta.vars) {
    ind <- c(ind, which(pheno.data$Condition==i))
  }
  pdata.short <- pheno.data[ind, ]
  exprs.short <- exprs[,which(colnames(exprs) %in% rownames(pdata.short))]
  exprs.short <- exprs.short[,order(match(colnames(exprs.short), rownames(pdata.short)))]
  
  # Create design matrix
  design = model.matrix(~factor(pdata.short$Condition, levels=meta.vars), data=pdata.short)
  colnames(design) <- c(meta.vars[1], paste(meta.vars[1], "vs", meta.vars[2], sep=""))
  # Fit with linear models
  fit <- lmFit(exprs.short, design)
  fit <- eBayes(fit)
  # Get all the genes with logFC, p-values, no filtering
  degs <- topTable(fit, coef=paste(meta.vars[1], "vs", meta.vars[2], sep=""), adjust.method="fdr", number=nrow(fit))

  # Merge degs with expression matrix
  exprs.degs <- merge(degs, exprs.short, by="row.names")
  if (noannotate) {
    rownames(exprs.degs) <- exprs.degs[, 1]
    exprs.degs <- exprs.degs[, -1]
  } else {
    colnames(exprs.degs)[1] <- "ENTREZID"
    # Add information about gene names
    EntrezID_Symbol<-select(org.Hs.eg.db, exprs.degs$ENTREZID, c("SYMBOL", "GENENAME"))
    exprs.degs <- cbind(EntrezID_Symbol, exprs.degs)
    exprs.degs <- exprs.degs[,-4]
  }
  
  return(exprs.degs)
}

filterDEGS <- function(degs, pval, fc, adj) {
  # Filter by p-values
  if (missing(adj)) {
    degs <- degs[degs$adj.P.Val < pval,]
  } else if (adj==FALSE) {
    degs <- degs[degs$P.Value < pval,]
  } else if (adj==TRUE) {
    degs <- degs[degs$adj.P.Val < pval,]
  }
  
  # Sort by logFC
  degs <- degs[order(abs(degs$logFC), decreasing = TRUE),]
  # Filter by logFC
  degs <- degs[abs(degs$logFC) > fc,]  
  return (degs)
}

getUniqueProbesets <- function(exprs, platform) {
  require(WGCNA)
  ## Get probeset to entrezid mapping
  probesetsID <- rownames(exprs)
  probesetsID_EntrezID<-select(get(paste(platform, ".db", sep="")), probesetsID, "ENTREZID")
  
  ## Replace probesetsIDs with gene IDs in expression data matrix
  
  # Exclude NA probesets
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$ENTREZID!="NA"),]
  # Exclude probesets mapped to different genes simultaneously
  n_occur <- data.frame(table(probesetsID_EntrezID$PROBEID))
  uniques <- n_occur[n_occur$Freq == 1,]$Var1
  probesetsID_EntrezID <- probesetsID_EntrezID[which(probesetsID_EntrezID$PROBEID %in% uniques),]
  # Filter expression matrix based on left probesets
  exprs <- exprs[which(rownames(exprs) %in% probesetsID_EntrezID$PROBEID),]
  
  # Select one probeset among the probesets mapped to the same gene based on maximum average value across the samples
  collapsed = collapseRows(exprs, probesetsID_EntrezID$ENTREZID, probesetsID_EntrezID$PROBEID, method="MaxMean")  
  exprs <- collapsed$datETcollapsed
  
  return(exprs)
}