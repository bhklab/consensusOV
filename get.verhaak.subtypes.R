get.verhaak.subtypes <- function(expression.matrix, entrez.ids) {
  load("/Users/greg/repos/consensusOV/genesets.verhaak.RData")
	
	# For ssGSEA scores for the new samples, use the intersecting genes
	genesets <- lapply(genesets, function(x) intersect(x, fData(eset)$gene) )
  
	expression.matrix <- exprs(eset)
  rownames(expression.matrix) <- fData(eset)$gene
  
  ## Get ssGSEA subtype scores
	gsva.out <- gsva(expression.matrix, genesets, method="ssgsea", tau=0.75, parallel.sz=4, mx.diff=FALSE, ssgsea.norm=FALSE)
  gsva.out <- t(gsva.out)
  
  gsva.out <- apply(gsva.out, 2, function(x) ( x - min(x) ) / ( max(x) - min(x) ))
  
  ## Classify each sample according to the max ssGSEA subtype score
  
  subclasses <- apply(gsva.out, 1, function(x) colnames(gsva.out)[which.max(x)])
  
  subclasses <- factor(subclasses, levels=names(genesets))
  ## Append a new column for Verhaak subtypes
  return(list(Verhaak.subtypes=subclasses, gsva.out=gsva.out))
}
