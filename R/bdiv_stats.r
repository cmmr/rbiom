#' Test beta diversity vs categorical or numeric metadata.
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
#'     bdiv_stats(biom, regr = "BMI", stat.by = "==Body Site")
#'     
#'     # The R code used to compute the stats is in $code.
#'     tbl <- bdiv_stats(biom, stat.by = "Sex")
#'     tbl$code

bdiv_stats <- function (
    biom, stat.by = NULL, regr = NULL, bdiv = "Bray-Curtis", 
    weighted = TRUE, tree = NULL, within = NULL, between = NULL, 
    test = ifelse(is.null(regr), "means", "trends"), model = "lm", 
    trans = ifelse(is.null(regr), "none", "rank"), 
    split.by = NULL, level = 0.95, p.adj = "fdr" ) {
  
  
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
  # Sanity checks.
  #________________________________________________________
  with(params, {
    
    validate_bdiv(max = Inf)
    
    validate_meta('stat.by',  cmp = TRUE)
    validate_meta('regr', null_ok = TRUE)
    validate_meta('split.by', cmp = TRUE, max = Inf, null_ok = TRUE)
    
    # Validates and appends to `within` and `between`.
    validate_meta_cmp(c('stat.by', 'split.by'))
    
  })
  
  
  
  
  params$md <- with(params, c(stat.by, regr, split.by))
  params$df <- do.call(bdiv_table, fun_params(bdiv_table, params))
  
  with(params, {
    if (length(bdiv) > 1) { split.by %<>% c('.bdiv')
    } else                { df %<>% rename_response(paste0('.', bdiv)) }
  })
  
  stats <- do.call(stats_table, fun_params(stats_table, params))
  
  
  
  set_cache_value(cache_file, stats)
  
  return (stats)
}
