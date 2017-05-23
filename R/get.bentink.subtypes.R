#' Get ovarian cancer subtypes as defined by Bentink et al., 2012
#' 
#' @param expression.matrix A matrix of gene expression values with rows as genes, columns as samples.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of \code{expression.matrix}
#' @return A list with first value \code{Bentink.subtypes} containing a factor of subtype names;
#' and second value \code{angio} containing the ouput of \code{genefu::ovcAngiogenic}
#' @examples
#' library(Biobase)
#' library(genefu)
#' data(GSE14764.eset)
#' expression.matrix <- exprs(GSE14764.eset)
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.bentink.subtypes(expression.matrix, entrez.ids)
#' @references Bentink et al. \emph{Angiogenic mRNA and microRNA gene expression signature predicts
#' a novel subtype of serous ovarian cancer.} PloS one (2012).
#' @export

# from genefu
.geneid.map <-
function(geneid1, data1, geneid2, data2, verbose=FALSE) {
  
  nn <- names(geneid1)
  geneid1 <- as.character(geneid1)
  names(geneid1) <- nn
  nn <- names(geneid2)
  geneid2 <- as.character(geneid2)
  names(geneid2) <- nn
  if(is.null(names(geneid1))) { names(geneid1) <- dimnames(data1)[[2]] }
  if(!missing(data2) && is.null(names(geneid2))) { names(geneid2) <- dimnames(data2)[[2]] }
  if(!missing(data1) && !missing(geneid1) && !missing(geneid2)) {
    ## remove probes without any measurements
    na.ix <- apply(data1, 2, function(x) { return(all(is.na(x))) })
    data1 <- data1[ , !na.ix, drop=FALSE]
    geneid1 <- geneid1[!na.ix]
  } else { stop("data1, geneid1 and geneid2 parameters are mandatory!") }
  if(!missing(data2)) {
    ## remove probes without any measurements
    na.ix <- apply(data2, 2, function(x) { return(all(is.na(x))) })
    data2 <- data2[ , !na.ix, drop=FALSE]
    geneid2 <- geneid2[!na.ix]
  } else { data2 <- NULL }
  
  gix1 <- !is.na(geneid1)
  gix2 <- !is.na(geneid2)
  
  geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
  if(length(geneid.common) == 0) {
    warning("no gene ids in common!")
    return(list("geneid1"=NA, "data1"=NA, "geneid2"=NA, "data2"=NA))
  }
  
  ## dataset1
  ## probes corresponding to common gene ids
  gg <- names(geneid1)[is.element(geneid1, geneid.common)]
  gid <- geneid1[is.element(geneid1, geneid.common)]
  ## duplicated gene ids
  gid.dupl <- unique(gid[duplicated(gid)])
  gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
  ## unique gene ids
  gid.uniq <- gid[!is.element(gid, gid.dupl)]
  gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
  ## data corresponding to unique gene ids
  datat <- data1[ ,gg.uniq,drop=FALSE]
  ## data for duplicated gene ids
  if(length(gid.dupl) > 0) {
    if(verbose) { message("\ndataset1 duplicates...") }
    ## compute the standard deviation with a penalization on the number of missing values
    ## this should avoid selecting the most variant probe with a lot of missing values
    pena <- apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
    pena <- log((nrow(data1) + 1) / (pena + 1)) + 1
    #pena <- 1
    sdr <- drop(apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
    mysd <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "sd"=sdr)
    mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
    mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
    datat <- cbind(datat, data1[ , mysd[ , "probe"], drop=FALSE])
  }
  data1 <- datat
  geneid1 <- geneid1[dimnames(data1)[[2]]]
    
  #dataset2
  if(is.null(data2)) {
    #keep arbitrarily the first occurence of each duplicated geneid
    geneid2 <- geneid2[!duplicated(geneid2) & is.element(geneid2, geneid.common)]
  }
  else {
    ## probes corresponding to common gene ids
    gg <- names(geneid2)[is.element(geneid2, geneid.common)]
    gid <- geneid2[is.element(geneid2, geneid.common)]
    ## duplicated gene ids
    gid.dupl <- unique(gid[duplicated(gid)])
    gg.dupl <- names(geneid2)[is.element(geneid2, gid.dupl)]
    ## unique gene ids
    gid.uniq <- gid[!is.element(gid, gid.dupl)]
    gg.uniq <- names(geneid2)[is.element(geneid2, gid.uniq)]
    ## data corresponding to unique gene ids
    datat <- data2[ ,gg.uniq,drop=FALSE]
    ## data for duplicated gene ids
    if(length(gid.dupl) > 0) {
      if(verbose) { message("\ndataset2 duplicates...") }
      ## compute the standard deviation with a penalization on the number of missing values
      ## this should avoid selecting the most variant probe with a lotof missing values
      pena <- apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
      pena <- log((nrow(data2) + 1) / (pena + 1)) + 1
      #pena <- 1
      sdr <- drop(apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
      mysd <- cbind("probe"=gg.dupl, "gid"=geneid2[gg.dupl], "sd"=sdr)
      mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
      mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
      datat <- cbind(datat, data2[ , mysd[ , "probe"], drop=FALSE])
    }
    data2 <- datat
    geneid2 <- geneid2[dimnames(data2)[[2]]]
  }
  
  #same order for the two datasets
  rix <- match(geneid2, geneid1)
  geneid1 <- geneid1[rix]
  data1 <- data1[ ,rix,drop=FALSE]
  return(list("geneid1"=geneid1, "data1"=data1, "geneid2"=geneid2, "data2"=data2))
}

# from genefu
.sig.score <-
function(x, data, annot, do.mapping=FALSE, mapping, size=0, cutoff=NA, signed=TRUE, verbose=FALSE) {
  
  if(missing(data) || missing(annot)) { stop("data and annot parameters must be specified") }
  x <- as.data.frame(x, stringsAsFactors=FALSE)
  if(nrow(x) == 0) { stop("empty gene list!"); }

  myprobe <- as.character(x[ ,"probe"])
  mygid <- as.character(x[ ,"EntrezGene.ID"])
  mycoef <- as.numeric(x[ ,"coefficient"])
  names(mycoef) <- names(mygid) <- names(myprobe) <- myprobe

  nix <- order(abs(mycoef), decreasing=TRUE, na.last=NA)
  myprobe <- myprobe[nix]
  mygid <- mygid[nix]
  mycoef <- mycoef[nix]
   
   if(do.mapping) { ## mapping is requested
    gid1 <- mygid
    gid2 <- as.character(annot[ ,"EntrezGene.ID"])
    names(gid2) <- dimnames(annot)[[1]]
    ## remove missing and duplicated geneids from the gene list
    rm.ix <- is.na(gid1) | duplicated(gid1)
    gid1 <- gid1[!rm.ix]
  
    rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
    if(is.na(rr$geneid1[1])) {
      #no gene ids in common
      res <- rep(NA, nrow(data))
      names(res) <- dimnames(data)[[1]]
      gf <- c("mapped"=0, "total"=nrow(x))
      if(verbose) { message(sprintf("probe candidates: 0/%i", nrow(x))) }
      return(list("score"=res, "mapping"=gf, "probe"=cbind("probe"=NA, "EntrezGene.ID"=NA, "new.probe"=NA)))
    }
    nix <- match(rr$geneid2, mygid)
    myprobe <- myprobe[nix]
    mygid <- mygid[nix]
    mycoef <- mycoef[nix]
    gid1 <- rr$geneid2
    if(is.null(names(gid1))) { stop("problem with annotations!") }
    gid2 <- rr$geneid1
    if(is.null(names(gid2))) { stop("problem with annotations!") }
    data <- rr$data1
  
    #change the names of probes in x and data
    names(mycoef) <- names(mygid) <- mygid <- names(myprobe) <- myprobe <- as.character(gid1)
    dimnames(data)[[2]] <- as.character(gid2)
  } else { ## no mapping
    nix <- is.element(myprobe, dimnames(data)[[2]])
    myprobe <- myprobe[nix]
    mygid <- mygid[nix]
    mycoef <- mycoef[nix]
    gid1 <- gid2 <- mygid
    data <- data[ ,myprobe,drop=FALSE]
  }
  if(length(myprobe) == 0) {
    if(verbose) { message(sprintf("probe candidates: 0/%i", size)) }
    tt <- rep(NA, nrow(data))
    names(tt) <- dimnames(data)[[1]]
    return(list("score"=tt, "mapping"=c("mapped"=0, "total"=nrow(x)), "probe"=cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))))
  }
  
  if(size == 0 || size > nrow(x)) { size <- length(myprobe) }
  nix <- 1:size
  myprobe <- myprobe[nix]
  mygid <- mygid[nix]
  mycoef <- mycoef[nix]
  gid1 <- gid1[nix]
  gid2 <- gid2[nix]
  if(!is.na(cutoff)) {
    nix <- abs(mycoef) > cutoff
    myprobe <- myprobe[nix]
    mygid <- mygid[nix]
    mycoef <- mycoef[nix]
    gid1 <- gid1[nix]
    gid2 <- gid2[nix]
   }
  probe.candp <- myprobe[mycoef >= 0]
  probe.candn <- myprobe[mycoef < 0]
  gf <- length(myprobe)

  gf <- c("mapped"=gf, "total"=nrow(x))
  if(verbose) { message(sprintf("probe candidates: %i/%i",gf[1], gf[2])) }

  nprobe <- c(probe.candp, probe.candn)
  myw <- c("p"=length(probe.candp) / length(nprobe), "n"=length(probe.candn) / length(nprobe))
  res <- rep(0, nrow(data))
  
  if(signed) {
    ## consider only the sign of the coefficients
    if(length(probe.candp) > 0) { res <- myw["p"] * (apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=sum, na.rm=TRUE) / apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=function(x) { return(sum(!is.na(x))) })) }
    if(length(probe.candn) > 0) { res <- res - myw["n"] * (apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=sum, na.rm=TRUE) / apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=function(x) { return(sum(!is.na(x))) })) }
  } else {
    ## consider the exact value of the coefficients
    if(length(probe.candp) > 0) { res <- myw["p"] * (apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=function(x, y) { nix <- is.na(x); return(sum(x * y, na.rm=TRUE) / sum(y[!nix])) }, y=abs(mycoef[probe.candp]))) }
    if(length(probe.candn) > 0) { res <- res - myw["n"] * (apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=function(x, y) { nix <- is.na(x); return(sum(x * y, na.rm=TRUE) / sum(y[!nix])) }, y=abs(mycoef[probe.candn]))) }
  }
  return(list("score"=res, "mapping"=gf, "probe"=cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))))
}

# from genefu 
.rescale <-
function(x, na.rm=FALSE, q=0) {
  if(q == 0) {
    ma <- max(x, na.rm=na.rm)
    mi <- min(x, na.rm=na.rm)
  } else {
    ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
    mi <- quantile(x, probs=q/2, na.rm=na.rm)
  }
  xx <- (x - mi) / (ma - mi)
  attributes(xx) <- list("names"=names(x), "q1"=mi,"q2"=ma)
  return(xx)
}

# ovcAngiogenic from genefu
.ovcAngiogenic <-
function(data, annot, hgs, gmap=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), do.mapping=FALSE, verbose=FALSE) {
    gmap <- match.arg(gmap)
    if(missing(hgs)) { hgs <- rep(TRUE, nrow(data)) }
    if(do.mapping) {
        if(!is.element(gmap, colnames(annot))) { stop("gmap is not a column of annot!") }
        if(verbose) { message("the most variant probe is selected for each gene") }
        sigt <- sigOvcAngiogenic[order(abs(sigOvcAngiogenic[ ,"weight"]), decreasing=FALSE), ,drop=FALSE]
        sigt <- sigt[!duplicated(sigt[ ,gmap]), ,drop=FALSE]
        gid2 <- sigt[ ,gmap]
        names(gid2) <- rownames(sigt)
        gid1 <- annot[ ,gmap]
        names(gid1) <- colnames(data)
        rr <- .geneid.map(geneid1=gid1, data1=data, geneid2=gid2)
        data <- rr$data1
        annot <- annot[colnames(data), ,drop=FALSE]
        sigt <- sigt[names(rr$geneid2), ,drop=FALSE]
        pold <- colnames(data)
        pold2 <- rownames(sigt)
        colnames(data) <- rownames(annot) <- rownames(sigt) <- paste("geneid", annot[ ,gmap], sep=".")
        mymapping <- c("mapped"=nrow(sigt), "total"=nrow(sigOvcAngiogenic))
        myprobe <- data.frame("probe"=pold, "gene.map"=annot[ ,gmap], "new.probe"=pold2)
    } else {
        gix <- intersect(rownames(sigOvcAngiogenic), colnames(data))
        if(length(gix) < 2) { stop("data do not contain enough gene from the ovcTCGA signature!") }
        data <- data[ ,gix,drop=FALSE]
        annot <- annot[gix, ,drop=FALSE]
        mymapping <- c("mapped"=length(gix), "total"=nrow(sigOvcAngiogenic))
        myprobe <- data.frame("probe"=gix, "gene.map"=annot[ ,gmap], "new.probe"=gix)
        sigt <- sigOvcAngiogenic[gix, ,drop=FALSE]
    }
    
    #data(modelOvcAngiogenic)
    ss <- .sig.score(x=data.frame("probe"=colnames(data), "EntrezGene.ID"=annot[ ,gmap], "coefficient"=sigt[ ,"weight"]), data=data, annot=annot, do.mapping=FALSE, signed=TRUE)$score
    ## rescale only with the high grade, late stage, serous (hgs) patients
    rr <- .rescale(ss[hgs], q=0.05, na.rm=TRUE)
    ## rescale the whole dataset
    pscore <- ((ss - attributes(rr)$q1) / (attributes(rr)$q2 - attributes(rr)$q1) - 0.5) * 2
    emclust.ts <- mclust::estep(modelName="E", data=pscore, parameters=modelOvcAngiogenic)
    dimnames(emclust.ts$z) <- list(names(pscore), c("Angiogenic.proba", "nonAngiogenic.proba"))
    class.ts <- mclust::map(emclust.ts$z, warn=FALSE)
    names(class.ts) <- names(pscore)
    sbt.ts <- class.ts
    sbt.ts[class.ts == 1] <- "Angiogenic"
    sbt.ts[class.ts == 2] <- "nonAngiogenic"
    sbts <- data.frame("subtype.score"=pscore, "subtype"=sbt.ts, emclust.ts$z)
    prisk <- as.numeric(sbts[ ,"subtype"] == "Angiogenic")
  names(prisk) <- names(pscore) <- rownames(data)
  return (list("score"=pscore, "risk"=prisk, "mapping"=mymapping, "probe"=myprobe, "subtype"=sbts))
}



get.bentink.subtypes <-
function(expression.matrix, entrez.ids) {
  ## Classify new samples
  entrez.ids <- as.character(entrez.ids)
  expression.matrix.t <- t(expression.matrix)
  annot <- as.data.frame(entrez.ids)
  colnames(annot) <- "entrezgene"
  rownames(annot) <- annot$entrezgene
  colnames(expression.matrix.t) <- rownames(annot)
  angio <- .ovcAngiogenic(data = expression.matrix.t, annot=annot, gmap="entrezgene", do.mapping = TRUE)
  Bentink.subtypes <- angio$subtype$subtype
  return(list(Bentink.subtypes=Bentink.subtypes, angio=angio))
}
