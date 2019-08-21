#' Get ovarian cancer subtypes as defined by Hao et al., 2017
#'
#' @details Hao et al., 2017 derived a gene signature to predict the tissue 
#' of origin of ovarian tumors as either fallopian tube (FT) or ovarian surface
#' epithelium (OSE).
#' 
#' The authors found that expression patterns of tissue-specific genes, prognostic
#' genes, and molecular markers support a dualistic tissue origin of ovarian cancer, 
#' from either FT or OSE.
#'
#' The subtype classifier considers 112 signature genes including 37 genes
#' upregulated in FT and 75 genes upregulated in OSE. A score is computed
#' that is designed to range from 0 to 1 for FT tumors, while OSE tumors have a 
#' score ranging from -1 to 0.
#'
#' @param expression.matrix A matrix of gene expression values with genes as
#' rows, samples as columns.
#' @param entrez.ids A vector of Entrez Gene IDs, corresponding to the rows of
#' \code{expression.matrix}.
#' @return A list with first value \code{tissue} containing a factor
#' of subtype names (tissue of origin); and second value \code{score} containing
#' the tissue-of-origin score.
#' @examples
#' library(Biobase)
#' data(GSE14764.eset)
#' expression.matrix <- exprs(GSE14764.eset)
#' entrez.ids <- as.character(fData(GSE14764.eset)$EntrezGene.ID)
#' get.hao.subtypes(expression.matrix, entrez.ids)
#' @references Hao et al. (2017) Integrated analysis reveals tubal- and
#' ovarian-originated serous ovarian cancer and predicts differential therapeutic
#' responses. Clinical Cancer Research, 23:7400-11.
#' @author Ludwig Geistlinger 
#' @importFrom utils read.delim
#' @export
get.hao.subtypes <- function(expression.matrix, entrez.ids) 
{
    # get the signature 
    sig.file <- system.file("extdata/FT_vs_OSE_signature.txt", 
                                package="consensusOV")

    sig <- read.delim(sig.file, as.is = TRUE)
    sig <- split(sig[,"ENTREZID"], sig[,"TISSUE"])
    
    # warn if certain signature genes are not present in the data
    rownames(expression.matrix) <- entrez.ids
    not.ft <- sum(!(sig$FT %in% rownames(expression.matrix)))
    if(not.ft) warning(paste(not.ft, "of", 
                                length(sig$FT), "FT genes not present"))

    not.ose <- sum(!(sig$OSE %in% rownames(expression.matrix)))    
    if(not.ose) warning(paste(not.ose, "of", 
                                length(sig$OSE), "OSE genes not present"))

    # scoring for one sample at a time
    .scoreSample <- function(evals)
    {   
        evals <- sort(evals)
        .getRank <- function(genes) match(genes, names(evals)) / length(evals)
        r <- lapply(sig, .getRank)
        nna <- sum(!is.na(unlist(r)))
        rs <-  sum(r[[1]], na.rm = TRUE) + sum(1 - r[[2]], na.rm = TRUE) 
        s <- 2 * rs / nna - 1 
        return(s)
    }   

    # apply to each sample & assign tissue of origin by thresholding
    scores <- apply(expression.matrix, 2, .scoreSample)
    tissues <- factor(ifelse(scores > 0, "FT", "OSE"))
    tscores <- list(tissues = tissues, scores = scores)

    return(tscores)
}
