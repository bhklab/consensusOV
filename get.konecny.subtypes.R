get.konecny.subtypes <- function(expression.matrix, entrez.ids) {
  load("/Users/greg/repos/consensusOV/centroids.konecny.RData")
  
  expression.matrix <- t(scale(t(expression.matrix)))
  
	## Classify using nearest centroid with Spearman's rho
  intersecting.entrez.ids <- intersect(entrez.ids, rownames(centroids.konecny))
  
  centroids.konecny[rownames(centroids.konecny) %in% intersecting.entrez.ids,]
  centroids.konecny <- centroids.konecny[as.character(intersecting.entrez.ids),]
  expression.matrix <- expression.matrix[entrez.ids %in% intersecting.entrez.ids,]
  
  expression.matrix <- as.data.frame(expression.matrix)
  spearman.cc.vals <- sapply(centroids.konecny, function(x) sapply(expression.matrix, function(y) cor(x, y , method="spearman")))
  
  subclasses <- apply(spearman.cc.vals, 1, function(x) as.factor(colnames(spearman.cc.vals)[which.max(x)]))
  
  subclasses <- factor(subclasses, levels=colnames(centroids.konecny))
  
  return(list(Konecny.subtypes=subclasses, spearman.cc.vals=spearman.cc.vals))
}
