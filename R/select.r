#' Reduce samples to a specific list
#' 
#' @param biom  A BIOM object, as returned from \link{read.biom}. Objects of
#'        class \code{matrix} or \code{simple_triplet_matrix} are also accepted.
#'        
#' @param samples  Sample names, indices, or a logical vector identifying
#'        the samples to keep. The latter two should be based on the order of
#'        sample names given by \code{colnames(biom$counts)}.
#'        
#' @param nTop  Selects this number of samples, taking the sample with the most
#'        observations first, then the sample with the second-most observations,
#'        etc. Ties will be randomly ordered. If \code{nTop} is higher than the 
#'        number of samples in the dataset, the entire dataset will be returned. 
#'        See note.
#'        
#' @param nRandom  Randomly selects this number of samples. If higher than the
#'        number of samples in the dataset, the entire dataset will be returned.
#'        See note.
#'        
#' @param seed  Random seed, used when selecting \code{nRandom} samples.
#'        
#' @param fast  Should subsetting the taxa table and phylogenetic tree be 
#'        skipped? These slow steps are often not necessary. (Default: FALSE)
#' 
#' Note: Generally, you will specify only one of the filters: \code{samples},
#' \code{nTop}, or \code{nRandom}. However, specifying multiple filters is
#' allowed; they will be applied in the order listed above.
#' 
#' @return An object of the same class as \code{biom}, subsetted as requested.
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
select <- function (biom, samples=NULL, nTop=NULL, nRandom=NULL, seed=0, fast=FALSE) {
  
  #--------------------------------------------------------------
  # Accept multiple types of input: BIOM, slam matrix, R matrix
  #--------------------------------------------------------------
  if (is(biom, 'BIOM'))                         { res <- biom$counts
  } else if (is(biom, 'simple_triplet_matrix')) { res <- biom
  } else if (is(biom, 'matrix')) { res <- slam::as.simple_triplet_matrix(biom)
  } else { stop("Invalid object type passed to rbiom::select()") }
  
  
  if (!is.null(samples)) {
    res <- try(res[,samples], silent=TRUE)
    if (is(res, "try-error"))
      stop(simpleError(sprintf("Cannot select samples: %s", res)))
  }
  
  if (!is.null(nTop)) {
    nTop <- min(ncol(res), nTop)
    set.seed(seed)
    ranking <- rank(slam::col_sums(res), ties.method="random")
    res     <- res[,ranking > ncol(res) - nTop]
  }
  
  if (!is.null(nRandom)) {
    nRandom <- min(ncol(res), nRandom)
    set.seed(seed)
    res <- res[,sort(sample(seq_len(ncol(res)), nRandom))]
  }
  
  
  #--------------------------------------------------------------
  # Drop taxa with zero observations
  #--------------------------------------------------------------
  if (isFALSE(fast))
    res <- res[slam::row_sums(res) > 0, ]
  
  
  #--------------------------------------------------------------
  # Input was a matrix, so return a matrix.
  #--------------------------------------------------------------
  if (is(biom, 'simple_triplet_matrix')) return (res)
  if (is(biom, 'matrix'))                return (as.matrix(res))
  
  
  #--------------------------------------------------------------
  # Alter the rest of the BIOM to match this new subset
  #--------------------------------------------------------------
  samples         <- colnames(res)
  taxa            <- rownames(res)
  
  biom$counts     <- biom$counts[taxa, samples]
  biom$metadata   <- biom$metadata[samples,,drop=FALSE]
  biom$info$shape <- dim(biom$counts)
  biom$info$nnz   <- length(biom$counts$v)
  
  
  # These steps are noncritical and can be skipped for speed
  if (isFALSE(fast)) {
    
    biom$taxonomy <- biom$taxonomy[taxa,,drop=FALSE]
    
    if (!is.null(biom$phylogeny))
      biom$phylogeny <- rbiom::subtree(biom$phylogeny, taxa)
    
    if (!is.null(biom$sequences))
      biom$sequences <- biom$sequences[taxa]
  }
  
  
  return (biom)
}
