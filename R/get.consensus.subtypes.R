#' Get consensusOV ovarian cancer subtypes
#'
#' @param expression.matrix A matrix of gene expression values with rows as
#' genes, columns as samples.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of
#' \code{expression.matrix}
#' @param concordant.tumors.only Logical. Should the classifier trained only on 
#' tumors that are concordantly classified by Helland, Konecny, and Verhaak?
#' Defaults to TRUE.
#' @param remove.using.cutoff Specify whether to classify NA for samples that do
#'  not meet a margin cutoff
#' @param percentage.dataset.removed If remove.using.cutoff is TRUE, then
#' classify this percentage of samples to NA based on margin values
#' @param .training.dataset ExpressionSet containing the training data. 
#' Defaults to the pooled dataset across selected MetaGxOvarian datasets.
#' @param .dataset.names.to.keep Names of MetaGxOvarian datasets to use for
#' training
#' @return A list with first value \code{consensusOV.subtypes} containing a
#' factor of subtype names; and second value \code{rf.probs} containing a matrix
#'  of subtype probabilities
#' @examples
#' library(Biobase)
#' data(GSE14764.eset)
#' expression.matrix <- exprs(GSE14764.eset)
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.consensus.subtypes(expression.matrix, entrez.ids)
#' @importFrom randomForest randomForest
#' @importFrom stats cor ecdf predict quantile sd
#' @importFrom utils combn
#' @export
get.consensus.subtypes <- function(expression.matrix, entrez.ids,
    concordant.tumors.only=TRUE, remove.using.cutoff=FALSE,
    percentage.dataset.removed=0.75,
    .training.dataset=consensus.training.dataset.full,
    .dataset.names.to.keep=names(esets.rescaled.classified.filteredgenes)) 
{

  ### Load training data
  message("Loading training data")
  expression.matrix <- t(scale(t(expression.matrix)))
  entrez.ids <- as.character(entrez.ids)

  # use full training dataset or only specified datasets?
  dataset.names.to.keep <- .dataset.names.to.keep
  cond1 <- length(dataset.names.to.keep) == length(esets.rescaled.classified.filteredgenes)
  cond2 <- all(dataset.names.to.keep == names(esets.rescaled.classified.filteredgenes))
  if(cond1 && cond2) 
  {
    # Using all datasets
    # use consensustraining.dataset.full
    # This file is produced from classificationAcrossDatasets.Rnw
    training.dataset <- .training.dataset
  } else{
    esets.training <- esets.rescaled.classified.filteredgenes[dataset.names.to.keep]
    training.dataset <- dataset.merging(esets.training, method = "intersect",
                                    standardization = "none")
  }

  # restrict training to concordant tumors?   
  if(concordant.tumors.only)
  {    
    subtype.correspondances <- data.frame(
        Konecny=c("C1_immL", "C2_diffL", "C3_profL", "C4_mescL"),
        Verhaak=c("IMR", "DIF", "PRO", "MES"),
        Helland=c("C2", "C4", "C5", "C1")
    )
    cases.to.keep <-
      match(training.dataset$Konecny.subtypes, subtype.correspondances$Konecny) ==
      match(training.dataset$Verhaak.subtypes, subtype.correspondances$Verhaak) &
      match(training.dataset$Verhaak.subtypes, subtype.correspondances$Verhaak) ==
      match(training.dataset$Helland.subtypes, subtype.correspondances$Helland)
    training.dataset <- training.dataset[,cases.to.keep]
  }  

  train.labels <- training.dataset$Verhaak.subtypes
  levels(train.labels) <- paste0(levels(train.labels), "_consensus")

  intersecting.entrez.ids <- as.character(
      intersect(fData(training.dataset)$EntrezGene.ID, entrez.ids)
  )

  message("Training Random Forest...")

  eids <- fData(training.dataset)$EntrezGene.ID
  ind <- match(intersecting.entrez.ids, eids)
  train.expression.matrix <- t(exprs(training.dataset)[ind,])

  comb.mat <- combn(seq_along(intersecting.entrez.ids), 2)
  train.pairwise.matrix <- apply(comb.mat, 2, function(pair) 
  	train.expression.matrix[,pair[1]] > train.expression.matrix[,pair[2]])
  train.pairwise.vals <- as.data.frame(train.pairwise.matrix)

  rf.model <- randomForest(x=train.pairwise.vals, y=train.labels)

  ind <- match(intersecting.entrez.ids, entrez.ids)
  test.expression.matrix <- t(expression.matrix[ind,])

  comb.mat <- combn(seq_along(intersecting.entrez.ids), 2)
  test.pairwise.matrix <- apply(comb.mat, 2, function(pair) 
    test.expression.matrix[,pair[1]] > test.expression.matrix[,pair[2]])
  test.pairwise.vals <- as.data.frame(test.pairwise.matrix)

  my.predictions <- predict(rf.model, newdata = test.pairwise.vals)
  my.predictions.probs <- predict(rf.model,
                                  newdata = test.pairwise.matrix,
                                  type = 'prob')

  if(remove.using.cutoff) {
   	subtr <- apply(my.predictions.probs, 1, function(row) sort(row)[3])
	my.predictions.margins <- rowMax(my.predictions.probs) - subtr    
	ecdf.evaluated <- ecdf(my.predictions.margins)(my.predictions.margins)
    na.ind <- ecdf.evaluated < percentage.dataset.removed
    my.predictions[na.ind] <- NA
  }

  names(my.predictions) <- NULL
  my.predictions <- factor(my.predictions, levels=c("IMR_consensus",
                                                    "DIF_consensus",
                                                    "PRO_consensus",
                                                    "MES_consensus"))

  return(list(consensusOV.subtypes=my.predictions, rf.probs=my.predictions.probs))
}
