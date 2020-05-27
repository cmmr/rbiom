#' Reduce samples to a specific list
#' 
#' @param biom  A BIOM object, as returned from \link{read.biom}.
#' @param samples  Sample names, indices, or a logical vector identifying
#'     the samples to keep. The latter two should be based on the order of
#'     sample names given by \code{colnames(biom$counts)}.
#' @param nTop  Selects this number of samples, taking the sample with the most
#'     observations first, then the sample with the second-most observations,
#'     etc. If \code{nTop} is higher than the number of samples in the dataset, 
#'     the entire dataset will be returned. See note.
#' @param nRandom  Randomly selects this number of samples. If higher than the
#'     number of samples in the dataset, the entire dataset will be returned.
#'     See note.
#' @param seed  Random seed, used when selecting \code{nRandom} samples.
#' 
#' Note: Generally, you will specify only one of the filters: \code{samples},
#' \code{nTop}, or \code{nRandom}. However, specifying multiple filters is
#' allowed; they will be applied in the order listed above.
#' 
#' @return A \code{BIOM} object.
#' @export
#' @seealso \code{\link{subset}}
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ex1 <- select(biom, c('HMP14', 'HMP22', 'HMP03'))
#'     ex2 <- select(biom, c(32, 11, 28, 16, 46, 5))
#'     ex3 <- select(biom, 1:50 %% 6 == 0)
#'     ex4 <- select(biom, nRandom = 10)
#'     ex5 <- select(biom, nTop = 5)
#'     ex6 <- select(biom, samples = 10:40, nTop = 20, nRandom = 10)
#'
select <- function (biom, samples=NULL, nTop=NULL, nRandom=NULL, seed=0) {
  
  res <- biom$counts
  
  if (!is.null(samples)) {
    res <- try(res[,samples], silent=TRUE)
    if (is(res, "try-error"))
      stop(simpleError(sprintf("Cannot select samples: %s", res)))
  }
  
  if (!is.null(nTop)) {
    nTop <- min(ncol(res), nTop)
    res <- res[,rank(slam::col_sums(res)) > ncol(res) - nTop]
  }
  
  if (!is.null(nRandom)) {
    nRandom <- min(ncol(res), nRandom)
    set.seed(seed)
    res <- res[,sort(sample(seq_len(ncol(res)), nRandom))]
  }
  
  samples <- colnames(res)
  taxa    <- names(which(slam::row_sums(res) > 0))
  
  biom$counts    <- biom$counts[taxa, samples]
  biom$taxonomy  <- biom$taxonomy[taxa   ,,drop=FALSE]
  biom$metadata  <- biom$metadata[samples,,drop=FALSE]
  if (!is.null(biom$phylogeny)) {
    biom$phylogeny <- rbiom::subtree(biom$phylogeny, taxa)
  }
  if (!is.null(biom$sequences)) {
    biom$sequences <- biom$sequences[taxa]
  }
  biom$info$shape <- dim(biom$counts)
  biom$info$nnz   <- length(biom$counts$v)
  
  return (biom)
}