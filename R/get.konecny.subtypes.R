#' Get ovarian cancer subtypes as defined by Konecny et al., 2014
#' 
#' @param expression.matrix A matrix of gene expression values with rows as genes, columns as samples.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of \code{expression.matrix}
#' @return A list with first value \code{Konecny.subtypes} containing a factor of subtype names;
#' and second value \code{spearman.cc.vals} containing the Spearman correlation values per subtype
#' @examples
#' data(GSE14764.eset)
#' expression.matrix <- GSE14764.eset
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.konecny.subtypes(expression.matrix, entrez.ids)
#' @references Konecny et al. \emph{Prognostic and therapeutic relevance of molecular subtypes in high-grade 
#' serous ovarian cancer.} Journal of the National Cancer Institute (2014).
get.konecny.subtypes <-
function(expression.matrix, entrez.ids) {
  load("/Users/greg/repos/consensusOV/data/centroids.konecny.RData")
  
  expression.matrix <- t(scale(t(expression.matrix)))
  entrez.ids <- as.character(entrez.ids)
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
