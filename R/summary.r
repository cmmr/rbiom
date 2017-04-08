#' Summarize the contents of a BIOM object
#' 
#' @name summary
#' @param object  A BIOM object, as returned from \link{read.biom}.
#' @param ...  Not used.
#' @return A distance matrix.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     summary(biom)
#'

summary.BIOM <- function (object, ...) {
  
  shortlist <- function (x, n) {
    res <- head(x, n)
    res <- if (length(x) > n) c(res, "...") else res
    res <- paste(collapse=", ", res)
  }
  
  sampleNames <- shortlist(colnames(object$counts),   5)
  taxaNames   <- shortlist(rownames(object$counts),   5)
  taxaRanks   <- shortlist(colnames(object$taxonomy), 9)
  metadata    <- shortlist(colnames(object$metadata), 5)
  tree        <- ifelse(is.null(object$phylogeny), "Absent", "Present")
  
  cat(paste(sep="\n", 
    sprintf("%7.0f Samples:  (%s)", ncol(object$counts),   sampleNames),
    sprintf("%7.0f Taxa:     (%s)", nrow(object$counts),   taxaNames),
    sprintf("%7.0f Ranks:    (%s)", ncol(object$taxonomy), taxaRanks),
    sprintf("%7.0f Metadata: (%s)", ncol(object$metadata), metadata),
    sprintf("        Tree:     %s", tree)
  ))
  
}

