#' Reduce samples to a specific list
#' 
#' @param biom   A BIOM object, as returned from [read_biom()]. Objects of
#'        class \code{matrix} or \code{simple_triplet_matrix} are also accepted.
#'        
#' @param samples   Sample names, indices, or a logical vector identifying
#'        the samples to keep. The latter two should be based on the order of
#'        sample names given by \code{colnames(biom$counts)}.
#'        
#' @param nTop   Selects this number of samples, taking the sample with the most
#'        observations first, then the sample with the second-most observations,
#'        etc. Ties will be randomly ordered. If \code{nTop} is higher than the 
#'        number of samples in the dataset, the entire dataset will be returned. 
#'        See note.
#'        
#' @param nRandom   Randomly selects this number of samples. If higher than the
#'        number of samples in the dataset, the entire dataset will be returned.
#'        See note.
#'        
#' @param seed   Random seed, used when selecting \code{nRandom} samples.
#'        
#' @param fast   Should subsetting the phylogenetic tree and sequences be 
#'        skipped? These slow steps are often not necessary. (Default: FALSE)
#' 
#' Note: Generally, you will specify only one of the filters: \code{samples},
#' \code{nTop}, or \code{nRandom}. However, specifying multiple filters is
#' allowed; they will be applied in the order listed above.
#' 
#' @return An object of the same class as \code{biom}, subsetted as requested.
#' 
#' @seealso [sample_subset()]
#' 
#' @export
#' @examples
#'   \dontrun{
#'     library(rbiom)
#'     
#'     ex1 <- sample_select(hmp50, c('HMP14', 'HMP22', 'HMP03'))
#'     ex2 <- sample_select(hmp50, c(32, 11, 28, 16, 46, 5))
#'     ex3 <- sample_select(hmp50, 1:50 %% 6 == 0)
#'     ex4 <- sample_select(hmp50, nRandom = 10)
#'     ex5 <- sample_select(hmp50, nTop = 5)
#'     ex6 <- sample_select(hmp50, samples = 10:40, nTop = 20, nRandom = 10)
#'  }

sample_select <- function (biom, samples=NULL, nTop=NULL, nRandom=NULL, seed=0, fast=FALSE) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("sample_select", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Accept multiple types of input: BIOM, slam matrix, R matrix
  #________________________________________________________
  if (is(biom, 'BIOM'))                         { res <- biom$counts
  } else if (is(biom, 'simple_triplet_matrix')) { res <- biom
  } else if (is(biom, 'matrix')) { res <- slam::as.simple_triplet_matrix(biom)
  } else { stop("Invalid object type passed to sample_select()") }
  
  
  if (!is_null(samples)) {
    res <- try(res[,samples], silent=TRUE)
    if (is(res, "try-error"))
      stop(sprintf("Cannot select samples: %s", res))
  }
  
  if (!is_null(nTop)) {
    nTop <- min(ncol(res), nTop)
    set.seed(seed)
    ranking <- rank(slam::col_sums(res), ties.method="random")
    res     <- res[,ranking > ncol(res) - nTop]
  }
  
  if (!is_null(nRandom)) {
    nRandom <- min(ncol(res), nRandom)
    set.seed(seed)
    res <- res[,sort(sample(seq_len(ncol(res)), nRandom))]
  }
  
  
  #________________________________________________________
  # Drop taxa with zero observations
  #________________________________________________________
  res <- res[slam::row_sums(res) > 0, ]
  
  
  #________________________________________________________
  # Input was a matrix, so return a matrix.
  #________________________________________________________
  if (is(biom, 'simple_triplet_matrix')) return (res)
  if (is(biom, 'matrix'))                return (as.matrix(res))
  
  
  #________________________________________________________
  # Alter the rest of the BIOM to match this new subset
  #________________________________________________________
  samples         <- colnames(res)
  taxa            <- rownames(res)
  
  biom$counts     <- biom$counts[taxa, samples]
  biom$metadata   <- biom$metadata[samples,,drop=FALSE]
  biom$info$shape <- dim(biom$counts)
  biom$info$nnz   <- length(biom$counts$v)
  biom$taxonomy   <- biom$taxonomy[taxa,,drop=FALSE]
  
  
  # These steps are noncritical and can be skipped for speed
  if (isFALSE(fast)) {
    
    if (!is_null(biom$phylogeny))
      biom$phylogeny <- tree_subset(biom$phylogeny, taxa)
    
    if (!is_null(biom$sequences))
      biom$sequences <- biom$sequences[taxa]
  }
  
  
  #________________________________________________________
  # Attach sample_select() call to provenance tracking
  #________________________________________________________
  cl      <- match.call()
  cl[[1]] <- as.name("sample_select")
  cl[[2]] <- as.name("biom")
  for (i in seq_along(cl)[-(1:2)]) {
    cl[i] <- list(eval.parent(cl[[i]]))
  }
  names(cl)[[2]] <- ""
  
  attr(biom, 'history') <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("biom <- %s", deparse1(cl)) ))
  
  
  set_cache_value(cache_file, biom)
  return (biom)
}
