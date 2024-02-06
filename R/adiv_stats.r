

#' Test alpha diversity differences for significance.
#' 
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
#'     biom <- rarefy(hmp50)
#'       
#'     adiv_stats(biom, stat.by = "Body Site")
#'       
#'     adiv_stats(biom, stat.by = "Sex", test = "pw_means")
#'     
#'     adiv_stats(biom, stat.by = "Body Site", split.by = "Sex")
#'       

adiv_stats <- function (
    biom, stat.by, adiv = "Shannon", split.by = NULL,
    trans = "none", test = "means", p.adj = "fdr" ) {
  
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment())
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  with(params, {
    validate_adiv(max = Inf)
    validate_meta('stat.by')
    validate_meta('split.by', null_ok = TRUE, max = Inf)
  })
  
  
  params$md <- with(params, c(stat.by, split.by))
  params$df <- do.call(adiv_table, fun_params(adiv_table, params))
  
  with(params, {
    if (length(adiv) > 1) { split.by %<>% c('.adiv')
    } else                { df %<>% rename_response(paste0('.', adiv)) }
  })
  
  stats <- do.call(stats_table, fun_params(stats_table, params))
  
  
  return (stats)
}
