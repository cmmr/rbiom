
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
  
  with_cache("as_percent", environment(), NULL, local({
  
  
    if (!is(biom, 'BIOM'))
      stop (simpleError('In as_percent(), biom must be a BIOM-class object.'))
    
    if (attr(biom, 'rarefaction') == 1)
      return (biom)
    
    divisor <- if (is_rarefied(biom)) { attr(biom, 'rarefaction', exact = TRUE)
    } else                            { sample_sums(biom)[biom$counts$j] }
    
    biom$counts$v <- biom$counts$v / divisor
    attr(biom, 'rarefaction') <- 1
    
    return (biom)
    
  }))
}
