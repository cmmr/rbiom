#' Run after manually editing a BIOM object's content.
#'
#' @param biom  The \code{BIOM} object to repair.
#' 
#' @param prune  Remove taxa and samples with zero observations. (Default: TRUE)
#'
#' @param fast  Should subsetting the phylogenetic tree and sequences be 
#'        skipped? These slow steps are often not necessary. (Default: FALSE)
#'        
#' @return A \code{BIOM} object.
#' @export


repair <- function (biom, prune=TRUE, fast=FALSE) {
  
  #________________________________________________________
  # Sanity Pre-checks
  #________________________________________________________
  if (!is(biom, "BIOM"))
      stop(simpleError("Invalid BIOM object."))
  
  
  #________________________________________________________
  # Find minimal set of sample/taxa names
  #________________________________________________________
  sn <- intersect(colnames(biom[['counts']]), rownames(biom[['metadata']]))
  tn <- intersect(rownames(biom[['counts']]), rownames(biom[['taxonomy']]))
  
  if (isFALSE(fast)) {
    
    if (!is_null(biom[['phylogeny']]))
      tn <- intersect(tn, tips(biom[['phylogeny']]))
    
    if (!is_null(biom[['sequences']]))
      tn <- intersect(tn, names(biom[['sequences']]))
  }
  
  
  #________________________________________________________
  # Prune taxa/samples with zero observations
  #________________________________________________________
  if (prune) {
    biom[['counts']] <- biom[['counts']][tn,sn]
    tn <- tn[slam::row_sums(biom[['counts']]) > 0]
    sn <- sn[slam::col_sums(biom[['counts']]) > 0]
  }
  
  
  #________________________________________________________
  # Apply minimal set of sample/taxa names
  #________________________________________________________
  biom[['counts']]   <- biom[['counts']][tn,sn]
  biom[['metadata']] <- biom[['metadata']][sn,,drop=F]
  biom[['taxonomy']] <- biom[['taxonomy']][tn,,drop=F]
  
  if (isFALSE(fast)) {
    
    if (!is_null(biom[['phylogeny']])) {
      if (length(tn) > 1) {
        biom[['phylogeny']] <- subtree(biom[['phylogeny']], tn)
      } else {
        biom[['phylogeny']] <- NULL
      }
    }
  
  if (!is_null(biom[['sequences']]))
    biom[['sequences']] <- biom[['sequences']][tn]
  }
  
  
  #________________________________________________________
  # Update `info` values
  #________________________________________________________
  biom[['info']][['shape']] <- c(length(tn), length(sn))
  
  if (all(biom[['counts']][['v']] %% 1 == 0)) {
    biom[['info']][['matrix_element_type']] <- "int"
  } else {
    biom[['info']][['matrix_element_type']] <- "float"
  }
  
  
  #________________________________________________________
  # Sanity Post-checks
  #________________________________________________________
  if (length(tn) == 0) warning("repair(): All taxa have been dropped.")
  if (length(sn) == 0) warning("repair(): All samples have been dropped.")
  
  
  #________________________________________________________
  # Attach repair() call to provenance tracking
  #________________________________________________________
  cl      <- match.call()
  cl[[1]] <- as.name("repair")
  cl[[2]] <- as.name("biom")
  for (i in seq_along(cl)[-(1:2)]) {
    cl[i] <- list(eval.parent(cl[[i]]))
  }
  names(cl)[[2]] <- ""
  attr(biom, 'history') %<>% c(paste("biom <-", deparse1(cl)))
  
  
  return (biom)
}
