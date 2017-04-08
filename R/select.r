#' Reduce samples to a specific list
#' 
#' @param biom  A BIOM object, as returned from \link{read.biom}.
#' @param samples  Sample names, indices, or a logical vector identifying
#'     the samples to keep. The latter two should be based on the order of
#'     sample names given by \code{colnames(biom$counts)}.
#' @return A \code{BIOM} object.
#' @export
#' @seealso \code{\link{subset}}
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     ex1 <- select(biom, c('HMP14', 'HMP22', 'HMP03'))
#'     ex2 <- select(biom, c(32, 11, 28, 16, 46, 5))
#'     ex3 <- select(biom, 1:50 %% 6 == 0)
#'
select <- function (biom, samples) {
  
  res <- try(biom$counts[,samples], silent=TRUE)
  if (is(res, "try-error"))
    stop(simpleError(sprintf("Cannot select samples: %s", res)))
  
  samples <- colnames(res)
  taxa    <- names(which(slam::row_sums(res) > 0))
  
  biom$counts    <- biom$counts[taxa, samples]
  biom$taxonomy  <- biom$taxonomy[taxa   ,,drop=FALSE]
  biom$metadata  <- biom$metadata[samples,,drop=FALSE]
  biom$phylogeny <- ape::drop.tip(
    phy = biom$phylogeny, 
    tip = setdiff(biom$phylogeny$tip.label, taxa))
  biom$info$shape <- dim(biom$counts)
  biom$info$nnz   <- length(biom$counts$v)
  
  return (biom)
}