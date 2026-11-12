#' Combines rare_corrplot and rare_stacked into a single figure.
#'
#' @inherit rare_corrplot
#' @inherit documentation_default
#'
#' @family rarefaction
#' @family visualization
#'
#' @param rline   Where to draw the rarefaction depth line on both plots.
#'        Set to `TRUE` to show an auto-selected rarefaction depth, `FALSE`
#'        to not show a line, or an integer for a custom position.
#'        Default: `TRUE`
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     rare_multiplot(hmp50, stat.by = "Body Site")
#'

rare_multiplot <- function (
    biom, adiv = "Shannon", layers = "tc", rline = TRUE,
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE,
    test = "none", fit = "log", at = NULL, level = 0.95, p.adj = "fdr",
    transform = "none", alt = "!=", mu = 0, caption = TRUE, check = FALSE,
    cpus = n_cpus(), ... ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('rare_multiplot', params)
  schema_check('rare_multiplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('rare_multiplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    #________________________________________________________
    # One rarefaction depth for both panels.
    #________________________________________________________
    if (isTRUE(params$rline))
      params$rline <- suggest_rarefy_depth(biom$counts)


    #________________________________________________________
    # Build the two subplots, then arrange them together.
    #________________________________________________________
    corr <- local({
      pp <- rlang::env_clone(params)
      pp$.fn <- 'rare_multiplot()'
      panel_params(pp, 1)
      rare_corrplot_build(pp)
    })

    stacked <- local({
      pp <- rlang::new_environment(parent = rlang::ns_env('rbiom'), data = list(
        biom        = biom,
        rline       = params$rline,
        counts      = TRUE,
        labels      = TRUE,
        y.transform = "log10",
        .dots       = list(),
        .fn         = 'rare_multiplot()' ))
      panel_params(pp, 2)
      rare_stacked_build(pp)
    })

    p <- plot_wrap(list(corr, stacked), ncol = 1)

    attr(p, 'stats') <- attr(corr, 'stats', exact = TRUE)

    set_cache_value(cache_file, p)
  }

  return (plot_prov(p, cmd, biom))
}
