#' Run after manually editing a BIOM object's content.
#'
#' @param biom  The \code{BIOM} object to repair.
#' @param prune  Remove taxa and samples with zero observations. (Default: TRUE)
#' @return A \code{BIOM} object.
#' @export


repair <- function (biom, prune=TRUE) {
  
  #--------------------------------------------------------------
  # Sanity Pre-checks
  #--------------------------------------------------------------
  if (!is(biom, "BIOM"))
      stop(simpleError("Invalid BIOM object."))
  
  
  #--------------------------------------------------------------
  # Find minimal set of sample/taxa names
  #--------------------------------------------------------------
  sn <- intersect(colnames(biom[['counts']]), rownames(biom[['metadata']]))
  tn <- intersect(rownames(biom[['counts']]), rownames(biom[['taxonomy']]))
  
  if (!is.null(biom[['phylogeny']]))
    tn <- intersect(tn, tips(biom[['phylogeny']]))
  
  if (!is.null(biom[['sequences']]))
    tn <- intersect(tn, names(biom[['sequences']]))
  
  
  #--------------------------------------------------------------
  # Prune taxa/samples with zero observations
  #--------------------------------------------------------------
  if (prune) {
    biom[['counts']] <- biom[['counts']][tn,sn]
    tn <- tn[slam::row_sums(biom[['counts']]) > 0]
    sn <- sn[slam::col_sums(biom[['counts']]) > 0]
  }
  
  
  #--------------------------------------------------------------
  # Apply minimal set of sample/taxa names
  #--------------------------------------------------------------
  biom[['counts']]   <- biom[['counts']][tn,sn]
  biom[['metadata']] <- biom[['metadata']][sn,,drop=F]
  biom[['taxonomy']] <- biom[['taxonomy']][tn,,drop=F]
  
  if (!is.null(biom[['phylogeny']])) {
    if (length(tn) > 1) {
      biom[['phylogeny']] <- subtree(biom[['phylogeny']], tn)
    } else {
      biom[['phylogeny']] <- NULL
    }
  }
  
  if (!is.null(biom[['sequences']]))
    biom[['sequences']] <- biom[['sequences']][tn]
  
  
  #--------------------------------------------------------------
  # Update `info` values
  #--------------------------------------------------------------
  biom[['info']][['shape']] <- c(length(tn), length(sn))
  
  if (all(biom[['counts']][['v']] %% 1 == 0)) {
    biom[['info']][['matrix_element_type']] <- "int"
  } else {
    biom[['info']][['matrix_element_type']] <- "float"
  }
  
  
  #--------------------------------------------------------------
  # Sanity Post-checks
  #--------------------------------------------------------------
  if (length(tn) == 0) warning("repair(): All taxa have been dropped.")
  if (length(sn) == 0) warning("repair(): All samples have been dropped.")
  
  
  #--------------------------------------------------------------
  # Attach repair() call to provenance tracking
  #--------------------------------------------------------------
  cl      <- match.call()
  cl[[1]] <- as.name("repair")
  cl[[2]] <- as.name("biom")
  for (i in seq_along(cl)[-(1:2)]) {
    cl[[i]] <- eval.parent(cl[[i]])
  }
  names(cl)[[2]] <- ""
  attr(biom, 'history') %<>% c(paste("biom <-", deparse1(cl)))
  
  
  return (biom)
}
