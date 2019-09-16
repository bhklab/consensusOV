#' Get ovarian cancer subtypes
#'
#' @param expression.dataset Either a matrix of gene expression values with rows
#'  as genes, columns as samples;
#' or a BioBase::ExpressionSet object from MetaGxOvarian.
#' If \code{expression.dataset} is a matrix, then \code{entrez.ids} must have
#' length equal to the number of rows of \code{expression.dataset}.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of
#' \code{expression.dataset}
#' @param method The subtyping method to use
#' @param ... Optional parameters to be passed to the low level function
#' @return A list with first value \code{Konecny.subtypes} containing a factor
#' of subtype names; and second value \code{spearman.cc.vals} containing the
#' Spearman correlation values per subtype
#' 
#' @examples
#' library(Biobase)
#' data(GSE14764.eset)
#' expression.matrix <- exprs(GSE14764.eset)
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.subtypes(expression.matrix, entrez.ids, method="Konecny")
#' 
#' @importFrom methods is
#' 
#' @export
get.subtypes <-
function(expression.dataset, entrez.ids=NULL,
         method=c("consensusOV", "Helland", "Verhaak", "Konecny", "Bentink"),
         ...) {
	method <- match.arg(method)
	if(is.null(entrez.ids) && !is(expression.dataset, "ExpressionSet")) {
		stop("If entrez.ids are not provided, expression.dataset must be of class 'ExpressionSet' from MetaGxOvarian")
	}
	if(!is.null(entrez.ids) && !is(expression.dataset, "matrix")) {
		stop("If entrez.ids are provided, expression.dataset must be of class 'matrix'")
	}
	if(!is.null(entrez.ids) && length(entrez.ids) != nrow(expression.dataset)) {
		stop("The length of entrez.ids must be equal to the number of rows of expression.dataset")
	}
	dataFromMetaGx <- is.null(entrez.ids)

	if(dataFromMetaGx) {
		expression.matrix <- exprs(expression.dataset)
		entrez.ids <- as.character(fData(expression.dataset)$EntrezGene.ID)
	} else {
		expression.matrix <- expression.dataset
	}

	if(method=="consensusOV") {
	  subtypes.output <- get.consensus.subtypes(expression.matrix, entrez.ids, ...)
	} else if(method == "Helland") {
	  subtypes.output <- get.helland.subtypes(expression.matrix, entrez.ids, ...)
	} else if(method == "Verhaak") {
	  subtypes.output <- get.verhaak.subtypes(expression.matrix, entrez.ids, ...)
	} else if(method == "Konecny") {
	  subtypes.output <- get.konecny.subtypes(expression.matrix, entrez.ids, ...)
	} else if(method == "Bentink") {
	  subtypes.output <- get.bentink.subtypes(expression.matrix, entrez.ids, ...)
	}

	return(subtypes.output)
}
