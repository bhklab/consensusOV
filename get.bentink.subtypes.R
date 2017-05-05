get.bentink.subtypes <- function(expression.matrix, entrez.ids) {
  ## Classify new samples
  expression.matrix.t <- t(expression.matrix)
  annot <- fData(eset)
  colnames(annot)[which(colnames(annot) == "EntrezGene.ID")] <- "entrezgene"
  angio <- genefu::ovcAngiogenic(data = expression.matrix.t, annot=annot, gmap="entrezgene", do.mapping = TRUE)
  Bentink.subtypes <- angio$subtype$subtype
  return(list(Bentink.subtypes=Bentink.subtypes, angio=angio))
}
