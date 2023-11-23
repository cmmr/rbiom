#' Split an rbiom object by metadata, apply function, and return results in a 
#' list.
#' 
#' Calls \code{plyr::dlply} internally. Consider setting 
#' \code{otu_tree(biom) <- NULL} to speed up creation of subseted rbiom objects.
#' 
#' @inherit bdply params
#'           
#' @param FUN   The function to execute on each \code{biom} subset. \code{FUN} 
#'        may return any object, all of which will be returned in a named list.
#'        
#' @return A list with the function outputs.
#' 
#' 
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     blply(hmp50, "Sex", n_samples)
#'     
#'     blply(hmp50, c("Body Site", "Sex"), function (b) {
#'       adm <- adiv_matrix(b)[,c("Shannon", "Simpson")]
#'       apply(adm, 2L, mean)
#'     })
#'     
#'     iters <- list(w = c(TRUE, FALSE), d = c("bray", "euclid"))
#'     blply(hmp50, "Sex", iters = iters, function (b, w, d) {
#'       r <- range(bdiv_distmat(biom = b, bdiv = d, weighted = w))
#'       round(data.frame(min = r[[1]], max = r[[2]]))
#'     })
#'     
#'     
#'
blply <- function (biom, vars, FUN, ..., iters = list(), prefix = FALSE) {
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  if (!is(FUN, 'function')) stop("Please provide a function to FUN.")
  
  
  dots  <- list(...)
  iters <- expand.grid(iters, stringsAsFactors = FALSE)
  
  
  
  #________________________________________________________
  # Simple cases where we're not faceting by metadata.
  #________________________________________________________
  if (is.null(vars)) {
  
    if (nrow(iters) == 0) {
      result <- list(do.call(FUN, c(list(biom), dots)))
      
    } else {
      result <- plyr::alply(iters, 1L, function (iter) {
        do.call(FUN, c(list(biom), dots, as.list(iter)))
      })
    }
    
  } else {
    
    try(silent = TRUE, validate_biom(clone = FALSE))
    if (!is(biom, 'rbiom'))
      stop("Can't apply metadata partitions to non-rbiom object.")
    
    
    data <- sample_metadata(biom)
    validate_meta('vars', col_type = 'cat', max = Inf, null_ok = TRUE)
    
    result <- plyr::dlply(data, ply_cols(vars), function (df) {
      
      sub_biom <- sample_select(biom, as.character(df[['.sample']]))
      
      if (nrow(iters) == 0)
        return (do.call(FUN, c(list(sub_biom), dots)))
      
      plyr::alply(iters, 1L, function (iter) {
        do.call(FUN, c(list(sub_biom), dots, as.list(iter))) })
    })
    
    
    #________________________________________________________
    # Un-nest lists created by dlply(adply()) call.
    #________________________________________________________
    if (nrow(iters) > 0) {
      
      split_labels <- attr(result, 'split_labels', exact = TRUE)
      result       <- do.call(c, result)
      
      if (is.null(split_labels)) {
        attr(result, 'split_labels') <- iters
        
      } else {
        attr(result, 'split_labels') <- plyr::adply(
          .data    = split_labels, 
          .margins = 1L, 
          .fun     = function (x) iters )
      }
    }
    
  }
  
  
  if (isTRUE(prefix) && nrow(iters) > 0)
    attr(result, 'split_labels') <- local({
      df <- attr(result, 'split_labels')
      i  <- seq_len(ncol(iters)) + length(vars)
      colnames(df)[i] <- paste0(".", colnames(iters))
      return (df)
    })
  
  
  return (result)
}

