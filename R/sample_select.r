#' Reduce samples to a specific list
#' 
#' @inherit documentation_default
#' 
#' @family metadata
#' @family samples
#'        
#' @param samples   Sample names, indices, or a logical vector identifying
#'        the samples to keep. The latter two should be based on the order of
#'        sample names given by \code{sample_names(biom)}.
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
#' Note: Generally, you will specify only one of the filters: \code{samples},
#' \code{nTop}, or \code{nRandom}. However, specifying multiple filters is
#' allowed; they will be applied in the order listed above.
#' 
#' @return A \code{rbiom} object, subsetted as requested.
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

sample_select <- function (biom, samples=NULL, nTop=NULL, nRandom=NULL, seed=0) {
  
  validate_biom(clone = FALSE)
  
  params  <- eval_envir(environment())
  history <- append_history('biom', params)
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Apply each filtering parameter in order.
  #________________________________________________________
  
  biom   <- params$biom
  counts <- otu_matrix(biom, sparse = TRUE)
  
  if (!is_null(params$samples)) {
    counts <- try(counts[,params$samples], silent=TRUE)
    if (is(counts, "try-error"))
      stop(sprintf("Cannot select samples: %s", counts))
  }
  
  if (!is_null(params$nTop)) {
    nTop <- min(ncol(counts), params$nTop)
    set.seed(params$seed)
    ranking <- rank(slam::col_sums(counts), ties.method="random")
    counts  <- counts[,ranking > ncol(counts) - nTop]
  }
  
  if (!is_null(params$nRandom)) {
    nRandom <- min(ncol(counts), params$nRandom)
    set.seed(params$seed)
    counts <- counts[,sort(sample(seq_len(ncol(counts)), nRandom))]
  }
  
  
  #________________________________________________________
  # Alter the rest of the rbiom object to match new subset
  #________________________________________________________
  otu_matrix(biom) <- counts
  
  
  
  attr(biom, 'display') <- "biom"
  attr(biom, 'history') <- history
  
  invalidate_biom()
  set_cache_value(cache_file, biom)
  
  return (biom)
}
