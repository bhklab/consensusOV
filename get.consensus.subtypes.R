
get.consensus.subtypes <- function(expression.matrix, entrez.ids, .dataset.names.to.keep=names(esets.not.rescaled.classified), remove.using.cutoff=FALSE, percentage.dataset.removed = 0.75) {
  
  ### Load training data
  print("Loading training data")
  ## This file is produced from classificationAcrossDatasets.Rnw
  expression.matrix <- t(scale(t(expression.matrix))) 
  
  load("/Users/greg/repos/consensusOV/entrez.id.union.RData")
  load("/Users/greg/repos/consensusOV/verhaak.entrez.ids.RData")
  
  # Save two datasets: the full merged dataset for classifying new cases, and a list of individual datasets for selecting datasets to use for training
  # 
  load("consensus.training.dataset.full.RData")
  load("esets.rescaled.classified.filteredgenes.RData")
  
  dataset.names.to.keep <- .dataset.names.to.keep
  
  if(length(dataset.names.to.keep) == length(esets.rescaled.classified.filteredgenes) &&
     all(dataset.names.to.keep == names(esets.rescaled.classified.filteredgenes))) {
    # Using all datasets
    # use consensustraining.datasetfull
    training.dataset <- consensus.training.dataset.full
  } else{
    esets.training <- esets.rescaled.classified.filteredgenes[dataset.names.to.keep]
    esets.merged <- MetaGx::datasetMerging(esets.training, method = "intersect", standardization = "none")
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
  
  rf.model <- randomForest(x=train.pairwise.vals, y=train.labels)
  
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
  
  return(list(consensusSubtpyes=my.predictions, rf.probs=my.predictions.probs))
}