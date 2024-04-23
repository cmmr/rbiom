#' Combines rare_corrplot and rare_stacked into a single figure.
#' 
#' @inherit rare_corrplot
#' @inherit documentation_default
#' 
#' @family rarefaction
#' @family visualization
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     rare_multiplot(hmp50, stat.by = "Body Site")
#'     

rare_multiplot <- function (
    biom, adiv = "Shannon", layers = "tc", rline = TRUE,
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = FALSE, 
    test = "none", fit = "log", at = NULL, level = 0.95, p.adj = "fdr", 
    trans = "none", alt = "!=", mu = 0, caption = TRUE, ... ) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('rare_multiplot', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Build the two subplots, then arrange them together.
  #________________________________________________________
  corrplot <- do.call(rare_corrplot, fun_params(rare_corrplot, params))
  stacked  <- do.call(rare_stacked,  fun_params(rare_stacked,  params))
  plots    <- list(corrplot = corrplot, stacked = stacked)
  
  
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
  
  
  
  attr(p, 'cmd') <-  current_cmd('rare_multiplot')
  set_cache_value(cache_file, p)
  
  return (p)
}






