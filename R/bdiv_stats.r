#' Test beta diversity vs categorical metadata.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_dist_test
#' @inherit documentation_default
#' @inherit documentation_stats_return return
#' 
#' @family beta_diversity
#' @family stats_tables
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     bdiv_stats(biom, stat.by = "Body Site", split.by = "Sex")
#'       
#'     bdiv_stats(biom, stat.by = "Sex", bdiv = c("bray", "unifrac"))
#'     
#'     # The R code used to compute the stats is in $code.
#'     tbl <- bdiv_stats(biom, stat.by = "Sex")
#'     tbl$code

bdiv_stats <- function (
    biom, stat.by, bdiv = "Bray-Curtis", weighted = TRUE, 
    split.by = NULL, within = NULL, between = NULL, tree = NULL, 
    trans = "none", test = "means", p.adj = "fdr" ) {
  
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment())
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  with(params, {
    
    validate_bdiv(max = Inf)
    
    validate_meta('stat.by',  cmp = TRUE)
    validate_meta('split.by', cmp = TRUE, max = Inf, null_ok = TRUE)
    
    # Validates and appends to `within` and `between`.
    validate_meta_cmp(c('stat.by', 'split.by'))
    
  })
  
  
  
  
  params$md <- with(params, c(stat.by, split.by))
  params$df <- do.call(bdiv_table, fun_params(bdiv_table, params))
  
  with(params, {
    if (length(bdiv) > 1) { split.by %<>% c('.bdiv')
    } else                { df %<>% rename_response(paste0('.', bdiv)) }
  })
  
  stats <- do.call(stats_table, fun_params(stats_table, params))
  
  
  
  set_cache_value(cache_file, stats)
  
  return (stats)
}
