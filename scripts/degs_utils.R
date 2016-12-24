getDEGS <- function(meta.vars, pheno.data, exprs) {

  ind <- c()
  for (i in meta.vars) {
    ind <- c(ind, which(pheno.data$Condition==i))
  }
  pdata.short <- pheno.data[ind, ]
  exprs.short <- exprs[, ind]
  #exprs.short <- exprs.short[,order(match(rownames(pdata.short),colnames(exprs)))]
  
  design = model.matrix(~as.factor(Condition), data=pdata.short)
  colnames(design) <- c(meta.vars[1], paste(meta.vars[1], "vs", meta.vars[2], sep=""))
  fit <- lmFit(exprs.short, design)
  fit <- eBayes(fit)
  degs <- topTable(fit, coef=paste(meta.vars[1], "vs", meta.vars[2], sep=""), adjust.method="fdr", number=nrow(fit))

  exprs.degs <- merge(degs, exprs.short, by="row.names")
  colnames(exprs.degs)[1] <- "ENTREZID"
  EntrezID_Symbol<-select(org.Hs.eg.db, exprs.degs$ENTREZID, c("SYMBOL", "GENENAME"))

  exprs.degs <- cbind(EntrezID_Symbol, exprs.degs)
  exprs.degs <- exprs.degs[,-4]

  return(exprs.degs)
  
}

filterDEGS <- function(degs, pval, fc) {
  degs <- degs[degs$adj.P.Val < pval,]
  degs <- degs[order(abs(degs$logFC)),]
  degs <- degs[abs(degs$logFC) > fc,]  
  return (degs)
}