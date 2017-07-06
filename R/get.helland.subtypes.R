#' Get ovarian cancer subtypes as defined by Helland et al., 2011
#'
#' @param expression.matrix A matrix of gene expression values with rows as
#' genes, columns as samples.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of
#' \code{expression.matrix}
#' @return A list with first value \code{Helland.subtypes} containing a factor
#' of subtype names; and second value \code{subtype.scores} containing a matrix
#' of subtype scores
#' @examples
#' library(Biobase)
#' data(GSE14764.eset)
#' expression.matrix <- exprs(GSE14764.eset)
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.helland.subtypes(expression.matrix, entrez.ids)
#' @references Helland et al. \emph{Deregulation of MYCN, LIN28B and LET7 in a
#' molecular subtype of aggressive high-grade serous ovarian cancers.}
#' PloS one (2011).
#' @export
get.helland.subtypes <-
function(expression.matrix, entrez.ids) {

  entrez.ids <- as.character(entrez.ids)
  eids <- lapply(entrez.id.logFC.list.helland, function (x) x$ENTREZID)
  helland.gene.set <- Reduce(function(x,y) union(x, y), eids)
  intersecting.entrez.ids <- intersect(helland.gene.set, entrez.ids)

  # Only keep genes present in both the supplementary and this eset
  entrez.id.logFC.list <- lapply(
      entrez.id.logFC.list.helland,
      function(x) x[x$ENTREZID %in% intersecting.entrez.ids, ])

  subtype.scores <- sapply(entrez.id.logFC.list, function(x) {
    ordered.expression.subset <- expression.matrix[match(x$ENTREZID, entrez.ids),]
	res <- colSums(ordered.expression.subset * x$logFC)
    return(res)
  })

  old.rownames <- rownames(subtype.scores)
  # Scale to mean=0, variance=1
  subtype.scores <- apply(subtype.scores, 2, scale)
  rownames(subtype.scores) <- old.rownames

  ind <- apply(subtype.scores, 1, which.max)
  subclasses <- colnames(subtype.scores)[ind]	
  subclasses <- factor(subclasses, levels=c("C2", "C4", "C5", "C1"))

  return(list(Helland.subtypes=subclasses, subtype.scores=subtype.scores))
}
