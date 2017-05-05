get.helland.subtypes <- function(expression.matrix, entrez.ids) {
  load("/Users/greg/repos/consensusOV/tothill.gene.set.RData")
  load("/Users/greg/repos/consensusOV/entrez.id.logFC.list.tothill.RData")
  
  intersecting.entrez.ids <- intersect(tothill.gene.set, entrez.ids)
  
  # Only keep genes present in both the supplementary and this eset
  entrez.id.logFC.list <- lapply(entrez.id.logFC.list.tothill, function(x) x[x$ENTREZID %in% intersecting.entrez.ids, ])
  
  subtype.scores <- sapply(entrez.id.logFC.list, function(x) {
    ordered.expression.subset <- exprs(eset)[match(x$ENTREZID, fData(eset)$EntrezGene.ID),]
    
    return(apply(ordered.expression.subset, 2, function(y) sum((y * x$logFC))))
    })
  old.rownames <- rownames(subtype.scores)
  # Scale to mean=0, variance=1
  subtype.scores <- apply(subtype.scores, 2, scale)
  rownames(subtype.scores) <- old.rownames
  subclasses <- factor(colnames(subtype.scores)[apply(subtype.scores, 1, which.max)], levels=names(supplementary.tables))
  
  return(list(Helland.subtypes=subclasses, subtype.scores=subtype.scores))
}
