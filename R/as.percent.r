
#' Convert absolute counts to relative abundances.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}.
#' @return A \code{BIOM} object.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     counts(hmp50)[1:4,1:5]
#'     
#'     biom <- as.percent(hmp50)
#'     counts(biom)[1:4,1:5]
#'

as.percent <- function (biom) {
  
  if (!is(biom, 'BIOM'))
    stop (simpleError('In as.percent(), biom must be a BIOM-class object.'))
  
  if (attr(biom, 'rarefaction') == 1)
    return (biom)
  
  divisor <- if (is.rarefied(biom)) { attr(biom, 'rarefaction', exact = TRUE)
  } else                            { sample.sums(biom)[biom$counts$j] }
  
  biom$counts$v <- biom$counts$v / divisor
  attr(biom, 'rarefaction') <- 1
  
  return (biom)
}
