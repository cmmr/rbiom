
#' Subset counts so that all samples have the same number of observations.
#'
#' @inherit documentation_default
#' 
#' @family rarefaction
#'     
#' @param depth   The number of observations to keep, per sample. If set to
#'        \code{NULL}, a depth will be automatically selected. Samples that 
#'        have fewer than this number of observations will be dropped. If 
#'        called on data with non-integer abundances, values will be re-scaled 
#'        to integers between 1 and \code{depth} such that they sum to 
#'        \code{depth}.
#'     
#' @param seed   An integer to use for seeding the random number generator. If
#'        you need to create different random rarefactions of the same 
#'        \code{rbiom} object, set this seed value to a different number each 
#'        time.
#'     
#' @return The rarefied \code{rbiom} object.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'
#'     range(sample_sums(hmp50))
#'
#'     biom <- sample_rarefy(hmp50, depth=1000)
#'     range(sample_sums(biom))
#'

sample_rarefy <- function (biom, depth=NULL, seed=0) {
  
  validate_biom(clone = FALSE)
  
  params  <- eval_envir(environment())
  history <- append_history('biom', params)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Choose the rarefaction depth to sample at
  #________________________________________________________
  if (is_null(depth)) depth <- rare_suggest(biom)
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  stopifnot(is_scalar_integerish(seed)  && !is_na(seed))
  stopifnot(is_scalar_integerish(depth) && !is_na(depth))
  stopifnot(depth > 0)
  
  

  counts <- otu_matrix(biom, sparse = TRUE)
  
  #________________________________________________________
  # Integer data. Randomly select observations to keep.
  #________________________________________________________
  
  if (all(counts$v %% 1 == 0)) {
    
    counts <- rcpp_rarefy(counts, depth, seed)
    
    
    #________________________________________________________
    # Fractional data. Re-scale between 1 and depth. Total = depth.
    #________________________________________________________
    
  } else {
    
    counts <- as.matrix(counts)
    
    # Does this, but then nudges up and down until == depth:
    # round(t((t(counts) / colSums(counts))) * depth)
    
    counts <- apply(counts, 2L, function (x) { 
      x        <- x / sum(x)
      y        <- round(x * depth)
      maxIters <- length(x)
      
      while (sum(y) < depth && maxIters > 0) {
        i        <- which.min( ((y + 1) / depth) - x )
        y[i]     <- y[i] + 1
        maxIters <- maxIters - 1
      }
      while (sum(y) > depth && maxIters > 0) {
        i        <- which.min( x - ((y - 1) / depth) )
        y[i]     <- y[i] - 1
        maxIters <- maxIters - 1
      }
      
      return (y)
    })
    
  }
  
  
  otu_matrix(biom) <- counts
  
  
  
  attr(biom, 'display')     <- "biom"
  attr(biom, 'rarefaction') <- depth
  attr(biom, 'history')     <- history
  
  invalidate_biom()
  set_cache_value(cache_file, biom)
  
  return (biom)
}




#' Suggest a 'good' rarefaction depth.
#'
#' @inherit documentation_default
#' 
#' @family rarefaction
#'     
#' @return An integer.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'
#'     rare_suggest(hmp50)
#'

rare_suggest <- function (biom) {
  
  validate_biom(clone = FALSE)
  
  
  if (all(biom$counts$v %% 1 == 0)) {
    
    sums  <- slam::col_sums(biom$counts)
    depth <- (sum(sums) * .1) / length(sums)
    depth <- min(sums[sums >= depth])
    
  } else {
    
    depth <- 10000
    
  }
  
  
  return (as.integer(depth))
}

