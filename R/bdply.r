#' Split BIOM by metadata, apply function, and return results in a data frame.
#' 
#' Calls \code{plyr::ddply} internally.
#' 
#' @name bdply
#' 
#' @param biom   A BIOM object, as returned from [read_biom()]. Technically
#'        accepts all types of objects, but only BIOM objects allow using the 
#'        \code{vars} option, below.
#' 
#' @param vars   A character vector of metadata fields. Each unique combination
#'        of values in these columns will be used to create a subsetted BIOM
#'        object to pass to \code{FUN}. If \code{NULL}, \code{biom} will be
#'        passed to \code{FUN} unaltered. Unambiguous abbreviations of metadata
#'        fields are also accepted.
#'           
#' @param FUN   The function to execute on each BIOM subset. \code{FUN} should
#'        return a data.frame, all of which will be \code{rbind}-ed together
#'        before being returned from \code{bdply}.
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
#' @param fast   If \code{TRUE} (the default), the subsetted BIOM objects will
#'        still contain the full taxa table and phylogenetic tree. Set 
#'        \code{fast = FALSE} to run the slow steps of subsetting these 
#'        elements as well.
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
#'       ad <- adiv_table(b, adiv = "all")[,c(".Shannon", ".Simpson")]
#'       apply(ad, 2L, mean)
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
bdply <- function (biom, vars, FUN, ..., iters = list(), prefix = FALSE, fast = TRUE) {
  
  
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
    
    if (!is(biom, 'BIOM'))
      stop("Can't apply metadata partitions to non-BIOM object.")
    
    data <- sample_metadata(biom, id = ".id")
    vars %<>% validate_arg(biom, 'meta', col_type = 'cat')
    
    result <- plyr::ddply(data, ply_cols(vars), function (df) {
      
      sub_biom <- sample_select(biom, df[['.id']], fast = fast)
      
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
  
  
  return (result)
  
}

