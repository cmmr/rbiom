#' Initialize caching for expensive biom operations.
#' 
#' Speeds up repetitive computations by storing results of functions calls. 
#' This cache initialization function can only be called once. Subsequent calls 
#' return a warning. 
#' 
#' Applies to \code{adiv_table()}, \code{apcoa()}, \code{bdply()}, 
#' \code{bdiv_dist()}, \code{bdiv_table()}, \code{distill()}, 
#' \code{ordinate()}, \code{plot()}, \code{rarefy()}, \code{select()}, 
#' \code{stats_table()}, \code{subtree()}, \code{taxa_rollup()}, and 
#' \code{unifrac()}.
#' 
#' @param cm  A memoise-compatible cache object. 
#'        E.g. \code{cachem::cache_mem()}, \code{cachem::cache_disk()}, or 
#'        \code{cachem::cache_layered()}. If \code{cm} is numeric, it will
#'        be used to define the maximum cache size (in bytes) for an in-memory 
#'        cache. Default: \code{50 * 1024^2} \emph{(50 MB)}
#' 
#' @param ...  Other arguments passed on to \code{cachem::cache_mem}. 
#'        For example: \code{max_size = 512 * 1024^2}, 
#'        \code{max_age = Inf}, \code{max_n = Inf}, etc.
#' 
#' @return The cache object, invisibly.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     init_cache()
#'     system.time(x <- adiv_table(hmp50, "multi"))
#'     system.time(y <- adiv_table(hmp50, "multi"))
#'     identical(x, y)
#'

init_cache <- function (cm = 50 * 1024^2, ...) {
  
  if (memoise::is.memoised(adiv_table))
    return(warning("rbiom is already memoised.\n"))
  
  if (is.numeric(cm)) cm <- cachem::cache_mem(max_size = cm, ...)
  if (is.null(cm))    cm <- cachem::cache_mem(...)
  
  assignInMyNamespace("adiv_table",  memoise::memoise(adiv_table,  cache = cm))
  assignInMyNamespace("apcoa",       memoise::memoise(apcoa,       cache = cm))
  assignInMyNamespace("bdply",       memoise::memoise(bdply,       cache = cm))
  assignInMyNamespace("bdiv_dist",   memoise::memoise(bdiv_dist,   cache = cm))
  assignInMyNamespace("bdiv_table",  memoise::memoise(bdiv_table,  cache = cm))
  assignInMyNamespace("distill",     memoise::memoise(distill,     cache = cm))
  assignInMyNamespace("ordinate",    memoise::memoise(ordinate,    cache = cm))
  assignInMyNamespace("plot.BIOM",   memoise::memoise(plot.BIOM,   cache = cm))
  assignInMyNamespace("rarefy",      memoise::memoise(rarefy,      cache = cm))
  assignInMyNamespace("select",      memoise::memoise(select,      cache = cm))
  assignInMyNamespace("stats_table", memoise::memoise(stats_table, cache = cm))
  assignInMyNamespace("subtree",     memoise::memoise(subtree,     cache = cm))
  assignInMyNamespace("taxa_rollup", memoise::memoise(taxa_rollup, cache = cm))
  assignInMyNamespace("unifrac",     memoise::memoise(unifrac,     cache = cm))
  
  return (invisible(cm))
}
