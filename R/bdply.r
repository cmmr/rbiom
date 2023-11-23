#' Split an rbiom object by metadata, apply function, and return results in a 
#' data frame.
#' 
#' Calls \code{plyr::ddply} internally. Consider setting 
#' \code{otu_tree(biom) <- NULL} to speed up creation of subseted rbiom objects.
#' 
#' @inherit documentation_default
#' 
#' @param vars   A character vector of metadata fields. Each unique combination
#'        of values in these columns will be used to create a subsetted 
#'        \code{rbiom} object to pass to \code{FUN}. If \code{NULL}, 
#'        \code{biom} will be passed to \code{FUN} unaltered. Unambiguous 
#'        abbreviations of metadata fields are also accepted.
#'           
#' @param FUN   The function to execute on each subset of \code{biom}. 
#'        \code{FUN} should return a data.frame, all of which will be 
#'        \code{rbind}-ed together before being returned from \code{bdply}.
#'                 
#' @param ...   Additional arguments to pass on to \code{FUN}.
#'           
#' @param iters   A named list of values to pass to \code{FUN}. Unlike 
#'        \code{...}, these will be iterated over in all combinations.
#'        Default: \code{list()}
#'           
#' @param prefix   When \code{TRUE}, prefixes the names in in \code{iters} with
#'        a '.' in the split_labels attribute of the returned object.
#'        Default: \code{FALSE}
#'        
#' @return A data.frame comprising the merged outputs of \code{FUN}, along with
#'         the columns specified by \code{vars}.
#' 
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     bdply(hmp50, "Sex", n_samples)
#'     
#'     bdply(hmp50, c("Body Site", "Sex"), function (b) {
#'       adm <- adiv_matrix(b)[,c("Shannon", "Simpson")]
#'       apply(adm, 2L, mean)
#'     })
#'     
#'     iters <- list(w = c(TRUE, FALSE), d = c("bray", "euclid"))
#'     bdply(hmp50, "Sex", iters = iters, function (b, w, d) {
#'       r <- range(bdiv_distmat(biom = b, bdiv = d, weighted = w))
#'       round(data.frame(min = r[[1]], max = r[[2]]))
#'     })
#'     
#'     
#'
bdply <- function (biom, vars, FUN, ..., iters = list(), prefix = FALSE) {
  
  
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
      result <- plyr::adply(iters, 1L, function (iter) {
        do.call(FUN, c(list(biom), dots, as.list(iter)))
      })
    }
    
  } else {
    
    
    
    try(silent = TRUE, validate_biom(clone = FALSE))
    if (!is(biom, 'rbiom'))
      stop("Can't apply metadata partitions to non-rbiom object.")
    
    data <- sample_metadata(biom)
    validate_meta('vars', col_type = 'cat', max = Inf, null_ok = TRUE)
    
    result <- plyr::ddply(data, ply_cols(vars), function (df) {
      
      sub_biom <- sample_select(biom, as.character(df[['.sample']]))
      
      if (nrow(iters) == 0)
        return (do.call(FUN, c(list(sub_biom), dots)))
      
      plyr::adply(iters, 1L, function (iter) {
        do.call(FUN, c(list(sub_biom), dots, as.list(iter))) })
    })
    
  }
  
  
  if (isTRUE(prefix) && nrow(iters) > 0) {
    i <- seq_len(ncol(iters)) + length(vars)
    colnames(result)[i] <- paste0(".", colnames(iters))
  }
  
  
  return (as_rbiom_tbl(result))
  
}

