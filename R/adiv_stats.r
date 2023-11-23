

#' Test alpha diversity differences for significance.
#' 
#' @inherit documentation_test.ifelse
#' @inherit documentation_model.lm
#' @inherit documentation_default
#' @inherit documentation_stats_return return
#' 
#' @family alpha_diversity
#' @family stats_tables
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- sample_rarefy(hmp50)
#'       
#'     adiv_stats(biom, stat.by = "Body Site")
#'       
#'     adiv_stats(biom, stat.by = "Sex", test = "pw_means")
#'       
#'     adiv_stats(biom, stat.by = "Body Site", regr = "Age")
#'     
#'     adiv_stats(biom, stat.by = "Body Site", split.by = "Sex")
#'       

adiv_stats <- function (
    biom, stat.by = NULL, regr = NULL, adiv = "Shannon",
    test = ifelse(is.null(regr), "means", "trends"), model = "lm", 
    split.by = NULL, level = 0.95, p.adj = "fdr") {
  
  
  validate_biom(clone = FALSE)
  
  params  <- eval_envir(environment())
  history <- append_history('stats ', params)
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  with(params, {
    
    validate_adiv(max = Inf)
    
    validate_meta('stat.by',  null_ok = TRUE)
    validate_meta('regr',     null_ok = TRUE)
    validate_meta('split.by', null_ok = TRUE, max = Inf)
    
  })
  
  
  params$md <- with(params, c(stat.by, regr, split.by))
  params$df <- do.call(adiv_table, fun_params(adiv_table, params))
  
  with(params, {
    if (length(adiv) > 1) { split.by %<>% c('.adiv')
    } else                { df %<>% rename_response(paste0('.', adiv)) }
  })
  
  stats <- do.call(stats_table, fun_params(stats_table, params))
  
  attr(stats, 'history') <- history
  return (stats)
}
