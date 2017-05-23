#' Get consensusOV ovarian cancer subtypes
#' 
#' @param expression.matrix A matrix of gene expression values with rows as genes, columns as samples.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of \code{expression.matrix}
#' @return A list with first value \code{consensusOV.subtypes} containing a factor of subtype names;
#' and second value \code{rf.probs} containing a matrix of subtype probabilities
#' @examples
#' library(Biobase)
#' data(GSE14764.eset)
#' expression.matrix <- exprs(GSE14764.eset)
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.consensus.subtypes(expression.matrix, entrez.ids)
#' @import randomForest
#' @export

# from MetaGx
.datasetMerging <- 
function (esets, method=c("union", "intersect"), standardization=c("quantile", "robust.scaling", "scaling", "none"), nthread=1) {
  
  method <- match.arg(method)
  standardization <- match.arg(standardization)
  
  ## all unique Entrez gene ids
  ## gene ids
  ugid <- lapply(esets, function(x) { return(Biobase::fData(x)) })
  ugid <- do.call(rbind, ugid)
  ugid <- ugid[!is.na(ugid[ , "EntrezGene.ID"]) & !duplicated(ugid[ , "EntrezGene.ID"]), , drop=FALSE]
  rownames(ugid) <- gsub(sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(ugid))
  switch (method,
    "union" = {
      feature.merged <- ugid
    },
    "intersect" = {
      feature.merged <- lapply(esets, function(x) { return(stripWhiteSpace(as.character(Biobase::fData(x)[ , "EntrezGene.ID"]))) })
      feature.merged <- table(unlist(feature.merged))
      feature.merged <- names(feature.merged)[feature.merged == length(esets)]
      feature.merged <- ugid[match(feature.merged, stripWhiteSpace(as.character(ugid[ , "EntrezGene.ID"]))), , drop=FALSE]
    },
    {
      stop("Unknown method")
    }
  )
  ## expression data
  exprs.merged <- lapply(esets, function (x, y) {
    ee <- Biobase::exprs(x)[is.element(rownames(exprs(x)),rownames(feature.merged)),]
    #print(dim(ee))
    eem <- matrix(NA, nrow=length(y), ncol=ncol(ee), dimnames=list(y, colnames(ee)))
    #print(dim(eem))
    #print(length(intersect(rownames(ee),rownames(eem))))
    eem[rownames(ee), colnames(ee)] <- ee
    return (eem)
  }, y=rownames(feature.merged))
  exprs.merged <- do.call(cbind, exprs.merged)
  ## clinical info
  ucid <- lapply(esets, function(x) { return(colnames(Biobase::pData(x))) })
  ucid <- table(unlist(ucid))
  ucid <- names(ucid)[ucid == length(esets)]
  clinicinfo.merged <- lapply(esets, function (x , y) {
    ee <- Biobase::pData(x)[ , y, drop=FALSE]
  }, y=ucid)
  clinicinfo.merged <- do.call(gdata::combine, clinicinfo.merged)
# if a data.source column already exists, remove it
  clinicinfo.merged$data.source <- NULL
  colnames(clinicinfo.merged)[which(colnames(clinicinfo.merged) == "source")] <- "data.source"
  rownames(clinicinfo.merged) <- colnames(exprs.merged)
#   rownames(clinicinfo.merged) <- gsub(   sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(clinicinfo.merged)         )
#   ## create a merged expressionSet object
  eset.merged <- ExpressionSet(assayData=exprs.merged, phenoData=AnnotatedDataFrame(data=clinicinfo.merged), featureData=AnnotatedDataFrame(data=feature.merged))
  experimentData(eset.merged)@preprocessing <- list("normalization"="mixed", package="unspecified", version="0")
  annotation(eset.merged) <- "mixed"
    
  ## standardization
  switch(standardization,
    "none" = {
      ## do nothing
    },
    {
      stop("Unknown data standardization method")
    }
  )
    
  return (eset.merged)
}


get.consensus.subtypes <-
function(expression.matrix, entrez.ids, .dataset.names.to.keep=names(esets.rescaled.classified.filteredgenes), remove.using.cutoff=FALSE, percentage.dataset.removed = 0.75) {
  
  ### Load training data
  print("Loading training data")
  ## This file is produced from classificationAcrossDatasets.Rnw
  expression.matrix <- t(scale(t(expression.matrix)))
  entrez.ids <- as.character(entrez.ids)
  
  
  dataset.names.to.keep <- .dataset.names.to.keep
  
  if(length(dataset.names.to.keep) == length(esets.rescaled.classified.filteredgenes) &&
     all(dataset.names.to.keep == names(esets.rescaled.classified.filteredgenes))) {
    # Using all datasets
    # use consensustraining.datasetfull
    training.dataset <- consensus.training.dataset.full
  } else{
    esets.training <- esets.rescaled.classified.filteredgenes[dataset.names.to.keep]
    esets.merged <- .datasetMerging(esets.training, method = "intersect", standardization = "none")
    subtype.correspondances <- data.frame(Konecny=c("C1_immL", "C2_diffL", "C3_profL", "C4_mescL"),
                                          Verhaak=c("IMR", "DIF", "PRO", "MES"),
                                          Helland=c("C2", "C4", "C5", "C1"))
    cases.to.keep <- 
      match(esets.merged$Konecny.subtypes, subtype.correspondances$Konecny) ==
      match(esets.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) &
      match(esets.merged$Verhaak.subtypes, subtype.correspondances$Verhaak) ==
      match(esets.merged$Helland.subtypes, subtype.correspondances$Helland)
    training.dataset <- esets.merged[,cases.to.keep]
  }
  
  ### Once we are happy with the normalization / removal of discordant cases, this eset should be a package data file.
  
  train.labels <- training.dataset$Verhaak.subtypes
  levels(train.labels) <- paste0(levels(train.labels), "_consensus")
  
  intersecting.entrez.ids <- as.character(intersect(fData(training.dataset)$EntrezGene.ID, entrez.ids))
  
  print("Training Random Forest...")
  
  train.expression.matrix <- t(exprs(training.dataset)[match(intersecting.entrez.ids, fData(training.dataset)$EntrezGene.ID),])
  
  train.pairwise.matrix <-
    apply(combn(1:length(intersecting.entrez.ids),2), 2, function(pair) train.expression.matrix[,pair[1]] > train.expression.matrix[,pair[2]])
  train.pairwise.vals <- as.data.frame(train.pairwise.matrix)
  
  rf.model <- randomForest::randomForest(x=train.pairwise.vals, y=train.labels)
  
  test.expression.matrix <- t(expression.matrix[match(intersecting.entrez.ids, entrez.ids),])
  
  test.pairwise.matrix <-
    apply(combn(1:length(intersecting.entrez.ids),2), 2, function(pair) test.expression.matrix[,pair[1]] > test.expression.matrix[,pair[2]])
  test.pairwise.vals <- as.data.frame(test.pairwise.matrix)
  
  my.predictions <- predict(rf.model, newdata = test.pairwise.vals)
  my.predictions.probs <- predict(rf.model, newdata = test.pairwise.matrix, type = 'prob')
  
  if(remove.using.cutoff) {
    my.predictions.margins <- rowMax(my.predictions.probs) - apply(my.predictions.probs, 1, function(row) sort(row)[3])
    my.predictions[ecdf(my.predictions.margins)(my.predictions.margins) < 0.75] <- NA
  }
  
  my.predictions <- factor(my.predictions, levels=c("IMR_consensus", "DIF_consensus", "PRO_consensus", "MES_consensus"))
  
  return(list(consensusOV.subtypes=my.predictions, rf.probs=my.predictions.probs))
}
