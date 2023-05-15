
#' Get the taxa abundances.
#' 
#' @name taxa_sums
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read_biom}.
#' @param rank  The taxonomic rank to return sums for. The default, 
#'        \code{NULL}, return OTU sums in the same order as they appear in 
#'        \link{counts}. If not \code{NULL}, returned sums are ordered from
#'        most abundance to least abundant.
#' @return A numeric vector, named with the taxa names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_sums(hmp50, 'Genus') %>% head(10)
#'

taxa_sums <- function (biom, rank = NULL) {
  
  with_cache(local({
  
    if (!is(biom, 'BIOM'))
      stop (simpleError('In taxa_sums(), biom must be a BIOM-class object.'))
    
    if (is_null(rank)) {
      x <- slam::row_sums(biom$counts)
    } else {
      x <- taxa_matrix(biom, rank) %>% colSums() %>% sort(decreasing = TRUE)
    }
    
    return (x)
    
  }))
}


#' Get the taxa abundances.
#' 
#' @name taxa_means
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read_biom}.
#' @param rank  The taxonomic rank to return means for. The default, 
#'        \code{NULL}, return OTU means in the same order as they appear in 
#'        \link{counts}. If not \code{NULL}, returned means are ordered from
#'        most abundance to least abundant.
#' @return A numeric vector, named with the taxa names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_means(hmp50, 'Genus') %>% head(10)
#'

taxa_means <- function (biom, rank = NULL) {
  
  with_cache(local({
    
    if (!is(biom, 'BIOM'))
      stop (simpleError('In taxa_means(), biom must be a BIOM-class object.'))
    
    if (is_null(rank)) {
      x <- slam::row_means(biom$counts)
    } else {
      x <- taxa_matrix(biom, rank) %>% colMeans() %>% sort(decreasing = TRUE)
    }
    
    return (x)
    
  }))
}


#' Get the names of the most abundant taxa.
#' 
#' @name top_taxa
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read_biom}.
#' @param rank  The taxonomic rank of interest. Default: \code{"OTU"}.
#' @param n  The number of taxa names to return. Default: \code{Inf}.
#' @return A character vector of the names of the top n most abundant taxa, 
#'         ordered from most abundant to least abundant.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     top_taxa(hmp50, 'Genus', 10)
#'

top_taxa <- function (biom, rank = 'OTU', n = Inf) {
  
  with_cache(local({
  
    if (!is(biom, 'BIOM'))
      stop (simpleError('In taxa_sums(), biom must be a BIOM-class object.'))
    
    x <- taxa_sums(biom, rank) %>% names() %>% head(n)
    
    return (x)
    
  }))
}
