
#' Run ordinations on a distance matrix.
#' 
#' @inherit documentation_default
#' 
#' @family ordination
#'        
#' @param ...  Additional arguments for \code{ord}.
#'        
#' @return A data.frame with columns \code{.sample}, \code{.ord}, \code{.x}, 
#'         \code{.y}, and (optionally) \code{.z}.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     dm  <- bdiv_distmat(hmp50, "bray")
#'     ord <- distmat_ord_table(dm, "PCoA")
#'     head(ord)
#'     

distmat_ord_table <- function (dm, ord = "PCoA", k = 2L, ...) {
  
  # browser() ## below code intermittently crashes R session
  
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env(..., .dots = TRUE)
  cache_file <- get_cache_file('distmat_ord_table', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  params <- list2env(params)
  with(params, {
    
    validate_ord(max = Inf)
    validate_var_range('k', range = c(2, 3), n = 1, int = TRUE)
    stopifnot(inherits(dm, 'dist'))
    
    if (is.null(attr(dm, 'display', exact = TRUE)))
      attr(dm, 'display') <- "dm"
    
  })
  
  
  
  
  #________________________________________________________
  # Call 3rd party function with proper arguments.
  #________________________________________________________
  dm   <- params$dm
  k    <- params$k
  ords <- setNames(params$ord, params$ord)
  dots <- params$.dots
  
  tbl <- plyr::ldply(ords, .id = ".ord", function (o) {
    
    if (o == "PCoA") {
      args <- c(fun_params(pcoa, dots), list(D = dm))
      res  <- do.call(pcoa, args)[['vectors']][,1:k]
      
    } else if (o == "tSNE") {
      require_package('tsne', 'to use "tSNE" ordination')
      args <- c(fun_params(tsne::tsne, dots), list(X = dm, k = k))
      res  <- suppressMessages(do.call(tsne::tsne, args))
      rownames(res) <- attr(dm, "Labels", exact = TRUE)
      
    } else if (o == "NMDS") {
      args <- c(fun_params(vegan::metaMDS, dots), list(comm = dm, k = k))
      args[['trace']] %<>% if.null(0)
      res  <- do.call(vegan::metaMDS, args)[['points']]
      rownames(res) <- attr(dm, "Labels", exact = TRUE)
      
    } else if (o == "UMAP") {
      require_package('uwot', 'to use "UMAP" ordination')
      args <- c(fun_params(uwot::umap, dots), list(X = dm))
      args[['n_components']] %<>% if.null(k)
      args[['n_neighbors']]  %<>% if.null(max(2, min(100, as.integer(attr(dm, 'Size') / 3))))
      res  <- do.call(uwot::umap, args)
    }
    
    params$.args <- args
    
    res <- as.data.frame(res)
    colnames(res) <- head(c(".x", ".y", ".z"), ncol(res))
    res[[".sample"]] <- rownames(res)
    rownames(res) <- NULL
    res %<>% keep_cols(".sample", ".x", ".y", ".z")
    
    return (res)
    
  }) %>% as_rbiom_tbl()
  
  
  
  
  #________________________________________________________
  # Table header.
  #________________________________________________________
  if (length(ords) == 1)
    attr(tbl, 'tbl_sum') <- c(
      'Ordination' = switch(
        EXPR = unname(ords),
        PCoA = "pcoa(%s)",
        tSNE = "tsne::tsne(%s)",
        NMDS = "vegan::metaMDS(%s)",
        UMAP = "uwot::umap(%s)" ) %>%
      sprintf(as.args(params$.args)) )
  
  
  attr(tbl, 'cmd') <- current_cmd('distmat_ord_table')
  
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}

