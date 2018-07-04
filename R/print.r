#' Summarize the contents of a BIOM object
#' 
#' @name print
#' @param x  A BIOM object, as returned from \link{read.biom}.
#' @param ...  Not used.
#' @return NULL (invisibly)
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     print(biom)
#'

print.BIOM <- function (x, ...) {
  
  shortlist <- function (object, n) {
    res <- head(object, n)
    res <- if (length(object) > n) c(res, "...") else res
    res <- paste(collapse=", ", res)
  }
  
  sampleNames <- shortlist(colnames(x$counts),   5)
  taxaNames   <- shortlist(rownames(x$counts),   5)
  taxaRanks   <- shortlist(colnames(x$taxonomy), 9)
  metadata    <- shortlist(colnames(x$metadata), 5)
  tree        <- ifelse(is.null(x$phylogeny), "Absent", "Present")
  
  cat(paste(sep="\n", 
    sprintf("%7.0f Samples:  (%s)", ncol(x$counts),   sampleNames),
    sprintf("%7.0f Taxa:     (%s)", nrow(x$counts),   taxaNames),
    sprintf("%7.0f Ranks:    (%s)", ncol(x$taxonomy), taxaRanks),
    sprintf("%7.0f Metadata: (%s)", ncol(x$metadata), metadata),
    sprintf("        Tree:     %s", tree), "\n"
  ))
  
  invisible(NULL)
}

