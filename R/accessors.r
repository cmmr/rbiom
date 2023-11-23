



#' Get \code{rbiom} object's miscellaneous information.
#' 
#' @family accessors
#' @inherit documentation_default
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
#' @family accessors
#' @inherit documentation_default
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
  return (sapply(as.list(biom[['metadata']][,-1]), is.numeric))
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
#' @family accessors
#' @inherit documentation_default
#' 
#' @return A named numeric vector of the number of observations in each 
#'         sample. The names are the sample IDs.
#'         
#' @seealso [rare_depth()], [rare_barplot()]
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
#' @return The rarefaction depth, or \code{NULL} if \code{biom} is not 
#'         rarefied.
#' 
#' @seealso [sample_sums()]
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


#' Check specific properties of an rbiom object.
#' 
#' @inherit documentation_default
#' 
#' @return \code{TRUE} if it has the property, \code{FALSE} otherwise.
#' @examples
#'     library(rbiom)
#'     
#'     has_tree(hmp50)
#'     has_sequences(hmp50)
#'     is_rarefied(hmp50)
#'     
#' @export
#' 
has_tree <- function (biom) {
  validate_biom(clone = FALSE)
  return (is(biom[['phylogeny']], 'phylo'))
}


#' @rdname has_tree
#' @export
has_sequences <- function (biom) {
  validate_biom(clone = FALSE)
  return (!is_null(biom[['sequences']]))
}


#' @rdname has_tree
#' @export
is_rarefied <- function (biom) {
  validate_biom(clone = FALSE)
  return (isTRUE(attr(biom, 'rarefaction', exact=TRUE) > 0))
}



#' Report the number of samples, otus, or ranks in an rbiom.
#' 
#' @family accessors
#' @inherit documentation_default
#' 
#' @return The number of samples, otus, taxa ranks, or metadata columns present.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     n_samples(hmp50)
#'     n_otus(hmp50)
#'     n_ranks(hmp50)
#'     n_metadata(hmp50)
#'

n_samples <- function (biom) {
  validate_biom(clone = FALSE)
  return (ncol(biom[['counts']]))
}


#' @rdname n_samples
#' @export

n_otus <- function (biom) {
  validate_biom(clone = FALSE)
  return (nrow(biom[['counts']]))
}


#' @rdname n_samples
#' @export

n_ranks <- function (biom) {
  validate_biom(clone = FALSE)
  n <- ncol(biom[['taxonomy']])
  return (n)
}


#' @rdname n_samples
#' @export

n_metadata <- function (biom) {
  validate_biom(clone = FALSE)
  n <- ncol(biom[['metadata']])
  return (n)
}
