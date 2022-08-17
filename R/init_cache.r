#' Initialize caching for expensive biom operations.
#' 
#' Speeds up repetitive computations by storing results of functions calls. 
#' This cache initialization function can only be called once. Subsequent calls 
#' return a warning.
#' 
#' @name init_cache
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
  
  fnList <- c(
    "adiv_boxplot", "adiv_table", "as_percent", "apcoa", "bdply", 
    "bdiv_boxplot", "bdiv_distmat", "bdiv_heatmap", "bdiv_biplot", 
    "bdiv_ord_table", "bdiv_table", "blply", "distill", "plot_heatmap", 
    "rarefy", "sample_sums", "select", "stats_table", "subtree", 
    "taxa_boxplot", "taxa_heatmap", "taxa_matrix", "taxa_means", 
    "taxa_stacked", "taxa_sums", "taxa_table", "top_taxa", "unifrac" )
  
  for (fn in fnList)
    assignInMyNamespace(fn, memoise::memoise(get(fn), cache = cm))
  
  return (invisible(cm))
}
