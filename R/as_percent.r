
#' Convert absolute counts to relative abundances.
#' 
#' @name as_percent
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read_biom}.
#' @return A \code{BIOM} object.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     counts(hmp50)[1:4,1:5]
#'     
#'     biom <- as_percent(hmp50)
#'     counts(biom)[1:4,1:5]
#'

as_percent <- function (biom) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("as_percent", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  if (!is(biom, 'BIOM'))
    stop (simpleError('In as_percent(), biom must be a BIOM-class object.'))
  
  if (attr(biom, 'rarefaction') == 1)
    return (biom)
  
  divisor <- if (is_rarefied(biom)) { attr(biom, 'rarefaction', exact = TRUE)
  } else                            { sample_sums(biom)[biom$counts$j] }
  
  biom$counts$v <- biom$counts$v / divisor
  attr(biom, 'rarefaction') <- 1
  
  set_cache_value(cache_file, biom)
  return (biom)
}
