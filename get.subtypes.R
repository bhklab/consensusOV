get.subtypes <- function(expression.dataset, entrez.ids=NULL, method=c("consensusOV", "Helland", "Verhaak", "Konecny"), ...) {
	method <- match.arg(method)
	if(is.null(entrez.ids) && class(expression.dataset != "ExpressionSet")) {
		stop("If entrez.ids are not provided, expression.dataset must be of class 'ExpressionSet' from MetaGxOvarian")
	}
	if(!is.null(entrez.ids) && class(expression.dataset != "matrix")) {
		stop("If entrez.ids are provided, expression.dataset must be of class 'matrix'")
	}
	if(!is.null(entrez.ids) && length(entrez.ids) != nrow(expression.dataset)) {
		stop("The length of entrez.ids must be equal to the number of rows of expression.dataset")
	}
	dataFromMetaGx <- is.null(entrez.ids)
	
	if(dataFromMetaGx) {
		expression.matrix <- exprs(expression.dataset)
		entrez.ids <- fData(esets.not.rescaled.classified$E.MTAB.386)$EntrezGene.ID
	} else {
		expression.matrix <- expression.dataset
	}

	subtypes.output
	
	if(method=="consensusOV") {
	  subtypes.output <- get.consensus.subtypes(expression.matrix, entrez.ids, ...)
	} else if(method == "Helland") {
	  subtypes.output <- get.helland.subtypes(expression.matrix, entrez.ids, ...)
	} else if(method == "Verhaak") {
	  subtypes.output <- get.verhaak.subtypes(expression.matrix, entrez.ids, ...)
	} else if(method == "Konecny") {
	  subtypes.output <- get.konecny.subtypes(expression.matrix, entrez.ids, ...)
	}
	
	return(subtypes.output)
}
