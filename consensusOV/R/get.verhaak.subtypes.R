#' Get ovarian cancer subtypes as defined by Verhaak et al., 2013
#' 
#' @param expression.matrix A matrix of gene expression values with rows as genes, columns as samples.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of \code{expression.matrix}
#' @return A list with first value \code{Verhaak.subtypes} containing a factor of subtype names;
#' and second value \code{gsva} containing the GSVA subtype scores
#' @examples
#' library(Biobase)
#' data(GSE14764.eset)
#' expression.matrix <- exprs(GSE14764.eset)
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.konecny.subtypes(expression.matrix, entrez.ids)
#' @references Verhaak et al. \emph{Prognostically relevant gene signatures of high-grade serous ovarian carcinoma.}
#' The Journal of Clinical Investigation (2013)
#' @importFrom GSVA gsva
#' @export
get.verhaak.subtypes <- function(expression.matrix, entrez.ids) {
  entrez.ids <- as.character(entrez.ids)
  rownames(expression.matrix) <- entrez.ids
  ## Get ssGSEA subtype scores
  gsva.out <- gsva(expression.matrix, verhaak.genesets.entrez.ids, method="ssgsea", tau=0.75, parallel.sz=4, mx.diff=FALSE, ssgsea.norm=FALSE)
  gsva.out <- t(gsva.out)
  
  gsva.out <- apply(gsva.out, 2, function(x) ( x - min(x) ) / ( max(x) - min(x) ))
  
  ## Classify each sample according to the max ssGSEA subtype score
  
  subclasses <- apply(gsva.out, 1, function(x) colnames(gsva.out)[which.max(x)])
  
  subclasses <- factor(subclasses, levels=c("IMR", "DIF", "PRO", "MES"))
  ## Append a new column for Verhaak subtypes
  return(list(Verhaak.subtypes=subclasses, gsva.out=gsva.out))
}
