#' Combines rare_corrplot and rare_barplot into a single figure.
#' 
#' @inherit documentation_corrplot
#' @inherit documentation_default
#' 
#' @family rarefaction
#' @family visualization
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     rare_multiplot(hmp50, color.by="Body Site")
#'     

rare_multiplot <- function (
    biom, adiv = "OTUs", layers = "t", rline = TRUE,
    color.by = NULL, facet.by = NULL, limit.by = NULL,
    formula = y ~ log(x), engine = "lm", level = 0.95, 
    caption = TRUE, labels = FALSE, ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("rare_multiplot(%s)", as.args(params, 2, rare_multiplot))
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Build the two subplots, then arrange them together.
  #________________________________________________________
  corrplot <- do.call(rare_corrplot, fun_params(rare_corrplot, params))
  barplot  <- do.call(rare_barplot,  fun_params(rare_barplot,  params))
  plots    <- list(corrplot = corrplot, barplot = barplot)
  
  
  p <- patchwork::wrap_plots(plots, ncol = 1) %>% 
    add_class('rbiom_plot')
  
  attr(p, 'stats') <- attr(corrplot, 'stats', exact = TRUE)
  attr(p, 'data')  <- lapply(plots, attr, which = 'data', exact = TRUE)
  
  attr(p, 'code') <- paste(collapse = "\n\n", local({
    
    cmds <- sapply(seq_along(plots), function (i) {
      sub(
        x           = attr(plots[[i]], 'code'), 
        pattern     = "ggplot(data", 
        replacement = sprintf("p%i <- ggplot(data[[%s]]", i, single_quote(names(plots)[[i]])),
        fixed       = TRUE )
    })
    c(cmds, sprintf("patchwork::wrap_plots(%s, ncol = 1)", paste0(collapse = ", ", "p", seq_along(plots))))
    
  })) %>% add_class('rbiom_code')
  
  
  
  attr(p, 'cmd') <- cmd
  set_cache_value(cache_file, p)
  
  return (p)
}






