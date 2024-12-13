#' Apply a function to each subset of an rbiom object.
#' 
#' `blply()` and `bdply()` let you divide your biom dataset into smaller 
#' pieces, run a function on those smaller rbiom objects, and return the 
#' results as a data.frame or list.
#' 
#' You can also specify additional variables for your function to iterate over
#' in unique combinations.
#' 
#' Calls [plyr::ddply()] or [plyr::dlply()] internally.
#' 
#' @inherit documentation_default
#' 
#' @family metadata
#' @family biom
#' 
#' @param vars   A character vector of metadata fields. Each unique combination
#'        of values in these columns will be used to create a subsetted 
#'        rbiom object to pass to `FUN.` If `NULL`, 
#'        `biom` will be passed to `FUN` unaltered. Unambiguous 
#'        abbreviations of metadata fields are also accepted.
#' 
#' @param FUN   The function to execute on each subset of `biom`. 
#'        For `bdply()`, the returned value will be coerced to a data.frame.
#'        For `blply()`, any returned value is unmodified.
#' 
#' @param ...   Additional arguments to pass on to `FUN`.
#'           
#' @param iters   A named list of values to pass to `FUN`. Unlike 
#'        \code{...}, these will be iterated over in all combinations.
#'        Default: \code{list()}
#' 
#' @param prefix   When `TRUE`, prefixes the names in in \code{iters} with
#'        a '.' in the final data.frame or 'split_labels' attribute.
#'        Default: `FALSE`
#' 
#' @return For `bdply()`, a tibble data.frame comprising the accumulated 
#'         outputs of `FUN`, along with the columns specified by 
#'         `vars` and `iters`. For `blply()`, a named list that has details
#'         about `vars` and `iters` in `attr(,'split_labels')`.
#' 
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     bdply(hmp50, "Sex", `$`, 'n_samples')
#'     
#'     blply(hmp50, "Sex", `$`, 'n_samples') %>% unlist()
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
bdply <- function (biom, vars, FUN, ..., iters = list(), prefix = FALSE) {
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  if (!inherits(FUN, 'function')) stop("Please provide a function to FUN.")
  
  
  dots  <- list(...)
  iters <- expand.grid(iters, stringsAsFactors = FALSE)
  
  
  
  #________________________________________________________
  # Simple cases where we're not faceting by metadata.
  #________________________________________________________
  if (is.null(vars)) {
    
    if (nrow(iters) == 0) {
      result <- do.call(FUN, c(list(biom), dots))
      
    } else {
      result <- plyr::adply(iters, 1L, function (iter) {
        do.call(FUN, c(list(biom), dots, as.list(iter)))
      })
    }
    
  } else {
    
    
    
    try(silent = TRUE, biom <- as_rbiom(biom))
    if (!inherits(biom, 'rbiom'))
      stop("Can't apply metadata partitions to non-rbiom object.")
    
    data <- biom$metadata
    validate_biom_field('vars', col_type = 'cat', max = Inf, null_ok = TRUE)
    
    result <- plyr::ddply(data, ply_cols(vars), function (df) {
      
      sub_biom        <- biom$clone()
      sub_biom$counts <- sub_biom$counts[,df[['.sample']]]
      
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



#' @rdname bdply
#' @export

blply <- function (biom, vars, FUN, ..., iters = list(), prefix = FALSE) {
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  if (!inherits(FUN, 'function')) stop("Please provide a function to FUN.")
  
  
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
    
    try(silent = TRUE, biom <- as_rbiom(biom))
    if (!inherits(biom, 'rbiom'))
      stop("Can't apply metadata partitions to non-rbiom object.")
    
    
    data <- biom$metadata
    validate_biom_field('vars', col_type = 'cat', max = Inf, null_ok = TRUE)
    
    result <- plyr::dlply(data, ply_cols(vars), function (df) {
      
      sub_biom        <- biom$clone()
      sub_biom$counts <- sub_biom$counts[,df[['.sample']]]
      
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

