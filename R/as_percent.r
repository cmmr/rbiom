
#' Convert absolute counts to relative abundances.
#' 
#' @inherit documentation_default
#' 
#' @family biom
#' 
#' @return An \code{rbiom} object.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     otu_matrix(hmp50)[1:4,1:5]
#'     
#'     biom <- as_percent(hmp50)
#'     otu_matrix(biom)[1:4,1:5]
#'

as_percent <- function (biom) {
  
  validate_biom(clone = TRUE)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment())
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  if (eq(attr(biom, 'rarefaction'), 1))
    return (biom)
  
  divisor <- if (is_rarefied(biom)) { attr(biom, 'rarefaction', exact = TRUE)
  } else                            { sample_sums(biom)[biom$counts$j] }
  
  biom$counts$v <- biom$counts$v / divisor
  attr(biom, 'rarefaction') <- 1
  
  
  invalidate_biom()
  set_cache_value(cache_file, biom)
  
  return (biom)
}
