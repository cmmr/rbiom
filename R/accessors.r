



#' Get \code{rbiom} object's miscellaneous information.
#' 
#' @inherit documentation_default
#' 
#' @family biom
#' 
#' @return A list of the top-level metadata in \code{biom}.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     biom_info(hmp50)
#'

biom_info <- function (biom) {
  validate_biom(clone = FALSE)
  return(biom[['info']])
}



#' Check which metadata columns are numeric.
#' 
#' @inherit documentation_default
#' 
#' @family metadata
#' 
#' @return A named logical vector indicating if each metadata column is 
#'         numeric.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     metadata_numeric(hmp50)
#'

metadata_numeric <- function (biom) {
  validate_biom(clone = FALSE)
  return (sapply(as.list(biom[['metadata']]), is.numeric))
}



#' Get the '.sample' metadata column's factor levels.
#' 
#' @noRd
#' @keywords internal
#' 
sample_levels <- function (biom) {
  validate_biom(clone = FALSE)
  return (levels(biom[['metadata']][['.sample']]))
}



#' Sum the observations in each sample.
#' 
#' @inherit documentation_default
#' 
#' @family samples
#' @family rarefaction
#' 
#' @return A named numeric vector of the number of observations in each 
#'         sample. The names are the sample IDs.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     sample_sums(hmp50) %>% sort() %>% head()
#'     
#'     hist(sample_sums(hmp50))

sample_sums <- function (biom) {
  validate_biom(clone = FALSE)
  return (slam::col_sums(biom[['counts']]))
}


#' If the samples are rarefied, report the level (depth).
#' 
#' This value is initally set by [read_biom()] and only subsequently changed by
#' calls to [sample_rarefy()]. Calling [as_percent()] will NOT set 
#' [rare_depth()] to \code{1}, even though that function sets all 
#' [sample_sums()] to \code{1}.
#' 
#' @inherit documentation_default
#' 
#' @family rarefaction
#' 
#' @return The rarefaction depth, or \code{NULL} if \code{biom} is not 
#'         rarefied.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     rare_depth(hmp50)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     rare_depth(biom)
#'

rare_depth <- function (biom) {
  validate_biom(clone = FALSE)
  return (attr(biom, 'rarefaction', exact = TRUE))
}


#' Check if a phylogenetic tree is included in an rbiom object.
#' 
#' @inherit documentation_default
#' 
#' @family phylogeny
#' 
#' @return \code{TRUE} if it has the property, \code{FALSE} otherwise.
#'     
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     has_tree(hmp50)
#' 
has_tree <- function (biom) {
  validate_biom(clone = FALSE)
  return (is(biom[['phylogeny']], 'phylo'))
}


#' Check if sequences are included in the rbiom object.
#' 
#' @inherit documentation_default
#' 
#' @family sequences
#' 
#' @return \code{TRUE} if \code{biom} has sequences, \code{FALSE} otherwise.
#'     
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     has_sequences(hmp50)
#' 
has_sequences <- function (biom) {
  validate_biom(clone = FALSE)
  return (!is_null(biom[['sequences']]))
}


#' Check if the rbiom object is rarefied.
#' 
#' @inherit documentation_default
#' 
#' @family rarefaction
#' 
#' @return \code{TRUE} if \code{biom} is rarefied, \code{FALSE} otherwise.
#'     
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     is_rarefied(hmp50)
#' 
is_rarefied <- function (biom) {
  validate_biom(clone = FALSE)
  return (isTRUE(attr(biom, 'rarefaction', exact=TRUE) > 0))
}



#' Report the number of samples in an rbiom object.
#' 
#' @inherit documentation_default
#' 
#' @family samples
#' 
#' @return The number of samples present.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     n_samples(hmp50)
#'
n_samples <- function (biom) {
  validate_biom(clone = FALSE)
  return (ncol(biom[['counts']]))
}


#' Report the number of otus in an rbiom object.
#' 
#' @inherit documentation_default
#' 
#' @family otus
#' 
#' @return The number of otus present.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     n_otus(hmp50)
#'     
n_otus <- function (biom) {
  validate_biom(clone = FALSE)
  return (nrow(biom[['counts']]))
}



#' Report the number of ranks in an rbiom object.
#' 
#' @inherit documentation_default
#' 
#' @family taxonomy
#' 
#' @return The number of taxonomic ranks present.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     n_ranks(hmp50)
#'
n_ranks <- function (biom) {
  validate_biom(clone = FALSE)
  n <- ncol(biom[['taxonomy']])
  return (n)
}



#' Report the number of metadata fields in an rbiom object.
#' 
#' @inherit documentation_default
#' 
#' @family metadata
#' 
#' @return The number of metadata fields present.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     n_metadata(hmp50)
#'
n_metadata <- function (biom) {
  validate_biom(clone = FALSE)
  n <- ncol(biom[['metadata']])
  return (n)
}
