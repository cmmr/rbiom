
#' Get summary taxa abundances.
#' 
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' 
#' @param rank  The taxonomic rank to return sums or means for. The default, 
#'        \code{0}, returns per-OTU summaries.
#'        
#' @return A named, sorted numeric vector.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     taxa_sums(hmp50) %>% head(4)
#'     
#'     taxa_means(hmp50, 'Family') %>% head(5)
#'

taxa_sums <- function (biom, rank = 0) {
  
  validate_biom(clone = FALSE)
  
  taxa_matrix(biom, rank, sparse = TRUE) %>%
    slam::row_sums() %>%
    sort(decreasing = TRUE)
}



#' @rdname taxa_sums
#' @export

taxa_means <- function (biom, rank = 0) {
  
  validate_biom(clone = FALSE)
  
  taxa_matrix(biom, rank, sparse = TRUE) %>%
    slam::row_means() %>%
    sort(decreasing = TRUE)
}
