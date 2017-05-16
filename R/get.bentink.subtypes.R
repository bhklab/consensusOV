#' Get ovarian cancer subtypes as defined by Bentink et al., 2012
#' 
#' @param expression.matrix A matrix of gene expression values with rows as genes, columns as samples.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of \code{expression.matrix}
#' @return A list with first value \code{Bentink.subtypes} containing a factor of subtype names;
#' and second value \code{angio} containing the ouput of \code{genefu::ovcAngiogenic}
#' @examples
#' library(Biobase)
#' data(GSE14764.eset)
#' expression.matrix <- GSE14764.eset
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.bentink.subtypes(expression.matrix, entrez.ids)
#' @references Bentink et al. \emph{Angiogenic mRNA and microRNA gene expression signature predicts
#' a novel subtype of serous ovarian cancer.} PloS one (2012).
#' @export
get.bentink.subtypes <-
function(expression.matrix, entrez.ids) {
  ## Classify new samples
  entrez.ids <- as.character(entrez.ids)
  expression.matrix.t <- t(expression.matrix)
  annot <- as.data.frame(entrez.ids)
  colnames(annot) <- "entrezgene"
  rownames(annot) <- annot$entrezgene
  colnames(expression.matrix.t) <- rownames(annot)
  angio <- genefu::ovcAngiogenic(data = expression.matrix.t, annot=annot, gmap="entrezgene", do.mapping = TRUE)
  Bentink.subtypes <- angio$subtype$subtype
  return(list(Bentink.subtypes=Bentink.subtypes, angio=angio))
}
