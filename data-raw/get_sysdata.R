library(Biobase)
library(MetaGxOvarian)
library(gdata)
library(hgu133plus2.db)
library(GSVA)
source("../R/dataset.merging.R")
source("../R/get.helland.subtypes.R")
source("../R/get.konecny.subtypes.R")
source("../R/get.verhaak.subtypes.R")
source("../R/get.bentink.subtypes.R")
source("../R/get.subtypes.R")

#### Get gene sets, coefficients from supplemental files ###

## centroids.konecny

konecny.supplementary.data <- gdata::read.xls(
    system.file("extdata", "jnci_JNCI_14_0249_s05.xls", package="MetaGx"),
    sheet=4, stringsAsFactors=FALSE)

# Subset supplementary.data to consist of centroids with intersecting genes
konecny.centroids <- konecny.supplementary.data[,c(2,4:7)]
konecny.centroids[,2:5] <- sapply(
    konecny.centroids[,2:5],
    function(x) ave(x, konecny.centroids$EntrezGeneID, FUN=mean))
konecny.centroids <- unique(konecny.centroids)
rownames(konecny.centroids) <- konecny.centroids$EntrezGeneID
konecny.centroids <- konecny.centroids[,-1]

## entrez.id.logFC.list.helland
helland.supplementary.type.1 <- gdata::read.xls(
    system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"),
    sheet=1, stringsAsFactors=FALSE)
helland.supplementary.type.2 <- gdata::read.xls(
    system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"),
    sheet=2, stringsAsFactors=FALSE)
helland.supplementary.type.4 <- gdata::read.xls(
    system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"),
    sheet=3, stringsAsFactors=FALSE)
helland.supplementary.type.5 <- gdata::read.xls(
    system.file("extdata", "journal.pone.0018064.s015.XLS", package="MetaGx"),
    sheet=4, stringsAsFactors=FALSE)
helland.supplementary.tables <- list(
    C1=helland.supplementary.type.1,
    C2=helland.supplementary.type.2,
    C4=helland.supplementary.type.4,
    C5=helland.supplementary.type.5)

entrez.id.logFC.list.helland <- lapply(helland.supplementary.tables,
                                       function(x) {
  ## Use the supplementary table's listed probe id and gene name to determine
  # the Entrez ID

  # If there is only one EntrezID that maps to a probe in hgu133plus2.db, use
  # that Entrez ID.

  # If there are multiple EntrezIDs that map to a probe, then use the EntrezID
  #(if any) that corresponds to the provided gene symbol.
  current.mapping <- suppressWarnings(
      AnnotationDbi::select(hgu133plus2.db, as.character(x$ID),
                            c("ENTREZID", "SYMBOL")))
  current.mapping <- current.mapping[ !is.na(current.mapping$ENTREZID), ]
  colnames(x)[1:2] <- c("PROBEID", "SYMBOL")
  mappings.with.unique.probeid <- current.mapping[ !(current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)]),]
  mappings.with.duplicate.probeid <- current.mapping[ current.mapping$PROBEID %in% current.mapping$PROBEID[duplicated(current.mapping$PROBEID)],]
  mappings.with.duplicate.probeid <- merge(x, mappings.with.duplicate.probeid, by=c("PROBEID", "SYMBOL"))[, c("PROBEID", "ENTREZID", "SYMBOL")]
  mappings.with.duplicate.probeid <- unique(mappings.with.duplicate.probeid)
  current.mapping <- rbind(mappings.with.unique.probeid, mappings.with.duplicate.probeid)
  to.return <- merge(x, current.mapping, by="PROBEID")[, c("ENTREZID", "PROBEID", "logFC")]
  return(to.return)
})

helland.gene.set <- Reduce(function(x,y) union(x, y),
                           lapply(entrez.id.logFC.list.helland,
                           function (x) x$ENTREZID))

## verhaak.genesets.entrez.ids
verhaak.supplementary.data.sheet7 <- read.xls(
    system.file("extdata", "JCI65833sd1.xls", package="MetaGx"),
    sheet=7, skip=1)

verhaak.gene.symbols <- lapply(levels(
    verhaak.supplementary.data.sheet7$CLASS),
    function(y) as.character(
        verhaak.supplementary.data.sheet7[
            verhaak.supplementary.data.sheet7$CLASS==y,1
        ]
    )
)


names(verhaak.gene.symbols) <- levels(verhaak.supplementary.data.sheet7$CLASS)

verhaak.genesets.entrez.ids <- lapply(
    verhaak.gene.symbols,
    function(gene.names) do.call(c, as.list(org.Hs.egALIAS2EG)[gene.names]))

## esets.rescaled.classified.filteredgenes
## consensus.training.dataset.full

verhaak.entrez.ids <- unlist(verhaak.genesets.entrez.ids)


##### Get classified datasets for consensus subtype classifier #####

remove.retracted <- FALSE
remove.subsets <- TRUE
probe.gene.mapping <- TRUE
keep.common.only <- FALSE
#only keep studies with at least this many samples
min.sample.size <- 40   
quantile.cutoff <- 0
strict.checking <- FALSE
remove.duplicates <- TRUE

rule.1 <- c("sample_type","^tumor$")
rule.2 <- c("histological_type","^ser$")
rule.3 <- c("summarystage","^late$")
rule.4 <- c("summarygrade","^high$")

rescale <- FALSE
# only keep esets with at least 10000 genes
min.number.of.genes <- 10000

source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))

esets.not.rescaled <- esets
# remove any genes with NA values
esets.not.rescaled <- lapply(esets.not.rescaled, function(eset) eset[apply(exprs(eset), 1, function(x) all(!is.na(x))), ])

esets.not.rescaled <- esets.not.rescaled[sapply(esets.not.rescaled, function(eset) nrow(eset) > 10000)]

# Classify data
set.seed(450)

# Subtype classification
if(file.exists("esets.not.rescaled.classified.RData")) {
  load("esets.not.rescaled.classified.RData")
} else {
  esets.not.rescaled.classified <- esets.not.rescaled
  # Remove TCGA.RNASeqV2
  esets.not.rescaled.classified <- esets.not.rescaled.classified[names(esets.not.rescaled.classified) != "TCGA.RNASeqV2"]
  
  esets.not.rescaled.classified <- lapply(esets.not.rescaled.classified, function(eset) {
    konecny.out <- get.subtypes(eset, method = "Konecny")
    annotated.eset <- eset
    annotated.eset$Konecny.subtypes <- konecny.out$Konecny.subtypes
    annotated.eset$Konecny.margins <- apply(konecny.out$spearman.cc.vals, 1, function(x) max(x) - sort(x)[3])
    return(annotated.eset)
  })
  esets.not.rescaled.classified <- lapply(esets.not.rescaled.classified, function(eset) {
    verhaak.out <- get.subtypes(eset, method = "Verhaak")
    annotated.eset <- eset
    annotated.eset$Verhaak.subtypes <- verhaak.out$Verhaak.subtypes
    annotated.eset$Verhaak.margins <- apply(verhaak.out$gsva.out, 1, function(x) max(x) - sort(x)[3])
    return(annotated.eset)
  })
  esets.not.rescaled.classified <- lapply(esets.not.rescaled.classified, function(eset) {
    helland.out <- get.subtypes(eset, method = "Helland")
    annotated.eset <- eset
    annotated.eset$Helland.subtypes <- helland.out$Helland.subtypes
    annotated.eset$Helland.margins <- apply(helland.out$subtype.scores, 1, function(x) max(x) - sort(x)[3])
    return(annotated.eset)
  })
  #esets.not.rescaled.classified <- lapply(esets.not.rescaled.classified, function(eset) getBentinkSubtypes(eset)[[1]])
  save(esets.not.rescaled.classified, file = "esets.not.rescaled.classified.RData")
}


esets.rescaled.classified <- lapply(esets.not.rescaled.classified,
                                    function(eset) {
                                        exprs(eset) <- t(scale(t(exprs(eset))))
                                        eset
                                    }
                            )

esets.rescaled.classified.filteredgenes <- lapply (
    esets.rescaled.classified,
    function(eset) {
        eset <- eset[
            as.character(fData(eset)$EntrezGene.ID) %in% as.character(verhaak.entrez.ids),
        ]
    }
)


consensus.training.dataset.full <- dataset.merging(
    esets.rescaled.classified.filteredgenes,
    method = "intersect",
    standardization = "none"
)

save(
  "konecny.centroids",
  "entrez.id.logFC.list.helland",
  "verhaak.genesets.entrez.ids",
  "esets.rescaled.classified.filteredgenes",
  "consensus.training.dataset.full",
  file="sysdata.rda", compress='xz')

## GSE14764.eset
GSE14764.eset <- esets.not.rescaled$GSE14764

# Only keep genes relevant for examples
entrez.id.union <- c(
    as.character(
        unlist(lapply(entrez.id.logFC.list.helland, function(x) x$ENTREZID))
    ),
    unlist(verhaak.genesets.entrez.ids),
    rownames(konecny.centroids)
)

GSE14764.eset <- GSE14764.eset[fData(GSE14764.eset)$EntrezGene.ID %in% entrez.id.union,]

save(GSE14764.eset, file="GSE14764.eset.RData", compress='xz')
