



#' Get \code{BIOM} object's miscellaneous information.
#' 
#' @param biom   A \code{BIOM} object, as returned from [read_biom()].
#' 
#' @return A list of the top-level metadata in \code{biom}.
#' 
#' @family accessors
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     biom_info(hmp50)
#'

biom_info <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (biom[['info']])
}



#' Sum the observations in each sample.
#' 
#' @family accessors
#' 
#' @param biom  A \code{BIOM} object, as returned from [read_biom()].
#' 
#' @return A named numeric vector of the number of observations in each 
#'         sample. The names are the sample IDs.
#'         
#' @seealso [rare_depth()]
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     sample_sums(hmp50) %>% sort() %>% head()
#'     
#'     hist(sample_sums(hmp50))

sample_sums <- function (biom, long=FALSE, md=FALSE) {
  
  stopifnot(is(biom, 'BIOM'))
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("sample_sums", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  result <- slam::col_sums(biom[['counts']])
  
  if (isTRUE(long) || !isFALSE(md)) {
    
    #________________________________________________________
    # Convert to long format
    #________________________________________________________
    result <- data.frame(
      stringsAsFactors = FALSE,
      'Sample' = names(result),
      'Reads'  = unname(result)
      )
    
    #________________________________________________________
    # Add Metadata
    #________________________________________________________
    if (identical(md, TRUE))  md <- colnames(sample_metadata(biom))
    if (identical(md, FALSE)) md <- c()
    for (i in unique(md))
      result[[i]] <- sample_metadata(biom, i)[result[['Sample']]]
    
  }
  
  set_cache_value(cache_file, result)
  return (result)
}


#' If the samples are rarefied, report the level (depth).
#' 
#' This value is initally set by [read_biom()] and only subsequently changed by
#' calls to [sample_rarefy()]. Calling [as_percent()] will NOT set 
#' [rare_depth()] to \code{1}, even though that function sets all 
#' [sample_sums()] to \code{1}.
#' 
#' @param biom  A \code{BIOM} object, as returned from [read_biom()].
#' 
#' @return The rarefaction depth, or \code{NULL} If the BIOM object is not 
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
  stopifnot(is(biom, 'BIOM'))
  return (attr(biom, 'rarefaction', exact = TRUE))
}


#' Check specific properties of a BIOM object.
#' 
#' @name properties
#' 
#' @param biom  A \code{BIOM} object, as returned from [read_biom()].
#' @return \code{TRUE} if it has the property, \code{FALSE} otherwise.
#' @examples
#'     library(rbiom)
#'     
#'     has_tree(hmp50)
#'     has_sequences(hmp50)
#'     is_rarefied(hmp50)
NULL


#' @rdname properties
#' @export
has_tree <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (is(biom[['phylogeny']], 'phylo'))
}


#' @rdname properties
#' @export
has_sequences <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (!is_null(biom[['sequences']]))
}


#' @rdname properties
#' @export
is_rarefied <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (isTRUE(attr(biom, 'rarefaction', exact=TRUE) > 0))
}



#' Report the number of samples, otus, or ranks in a BIOM.
#' 
#' @param biom   A \code{BIOM} object, as returned from [read_biom()].
#' 
#' @return The number of samples, otus, or ranks present.
#' 
#' @family accessors
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     n_samples(hmp50)
#'     n_otus(hmp50)
#'     n_ranks(hmp50)
#'

n_samples <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (ncol(biom[['counts']]))
}


#' @rdname n_samples
#' @family accessors
#' @export

n_otus <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (nrow(biom[['counts']]))
}


#' @rdname n_samples
#' @family accessors
#' @export

n_ranks <- function (biom) {
  stopifnot(is(biom, 'BIOM'))
  return (ncol(biom[['taxonomy']]))
}
