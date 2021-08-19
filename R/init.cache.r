#' Initialize caching for expensive biom operations.
#' 
#' Speeds up repetitive computations by storing results of functions calls. 
#' This cache initialization function can only be called once. Subsequent calls 
#' return an error. 
#' 
#' Applies to \code{alpha.div()}, \code{apcoa()}, \code{bdply()}, 
#' \code{beta.div()}, \code{distill()}, \code{ordinate()}, \code{plot()}, 
#' \code{rarefy()}, \code{select()}, \code{stats.table()}, \code{subtree()}, 
#' \code{taxa.rollup()}, and \code{unifrac()}.
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
#'     init.cache()
#'     system.time(x <- alpha.div(hmp50, "multi"))
#'     system.time(y <- alpha.div(hmp50, "multi"))
#'     identical(x, y)
#'

init.cache <- function (cm = 50 * 1024^2, ...) {
  
  if (memoise::is.memoised(alpha.div))
    stop("rbiom is already memoised.\n")
  
  if (is.numeric(cm)) cm <- cachem::cache_mem(max_size = cm, ...)
  if (is.null(cm))    cm <- cachem::cache_mem(...)
  
  assignInMyNamespace("alpha.div",   memoise::memoise(alpha.div,   cache = cm))
  assignInMyNamespace("apcoa",       memoise::memoise(apcoa,       cache = cm))
  assignInMyNamespace("bdply",       memoise::memoise(bdply,       cache = cm))
  assignInMyNamespace("beta.div",    memoise::memoise(beta.div,    cache = cm))
  assignInMyNamespace("distill",     memoise::memoise(distill,     cache = cm))
  assignInMyNamespace("ordinate",    memoise::memoise(ordinate,    cache = cm))
  assignInMyNamespace("plot.BIOM",   memoise::memoise(plot.BIOM,   cache = cm))
  assignInMyNamespace("rarefy",      memoise::memoise(rarefy,      cache = cm))
  assignInMyNamespace("select",      memoise::memoise(select,      cache = cm))
  assignInMyNamespace("stats.table", memoise::memoise(stats.table, cache = cm))
  assignInMyNamespace("subtree",     memoise::memoise(subtree,     cache = cm))
  assignInMyNamespace("taxa.rollup", memoise::memoise(taxa.rollup, cache = cm))
  assignInMyNamespace("unifrac",     memoise::memoise(unifrac,     cache = cm))
  
  return (invisible(cm))
}
