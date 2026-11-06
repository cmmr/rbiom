#________________________________________________________
# Computes p-values for regression plots and places them on the plot.
#
# Inputs read from `params`:
#   .ggdata, .xcol, .ycol, stat.by, facet.by, test, fit, at, level, alt,
#   mu, p.adj, p.top, caption, strata (optional), and, when
#   .ycol == ".distance" (beta diversity corrplots): biom, bdiv, tree,
#   alpha, cpus, seed, permutations.
#
# Which test runs:
#   .ycol == ".distance": a Mantel test (bdiv_stats(test = "mantel")) of
#       the distances against |x_i - x_j|, one per facet and stat.by
#       group. One label per facet; no vertical line.
#   otherwise: stats_table() with the plot's regr/resp/stat.by/facet.by,
#       honoring test ("emmeans"/"emtrends"), fit, at, level, alt, mu,
#       p.adj, and strata. Tests evaluated at an `at` position (pairwise
#       means; gam slopes) get a label and a dotted vertical line at that
#       position (the most significant row per facet when several); global
#       tests get a label in the top-left corner.
#
# Outputs:
#   params$.plot_attrs$stats     - the statistics table.
#   attr(.ggdata, 'stat_labels') - tibble: facet.by columns, .xcol, .ycol,
#                                  .label, .hjust, .vjust ('stats_label').
#   attr(.ggdata, 'stat_vline')  - tibble: facet.by columns, .xcol
#                                  ('stats_vline'); absent for global tests.
#________________________________________________________
corrplot_stats <- function (params) {

  stopifnot(is_bare_environment(params))


  # Validation of user arguments handled by corrplot_build()


  #________________________________________________________
  # Omit stats from this plot.
  #________________________________________________________
  if (eq(params$test, 'none'))
    return (params)

  is_bdiv <- eq(params$.ycol, '.distance')


  #________________________________________________________
  # If anything goes wrong, skip stats. (Optional parameters such as
  # `strata` or `biom` are NULL when the plot function doesn't set them.)
  #________________________________________________________
  stats <- tryCatch(
    error = function (e) { cli_warn(c('!' = "Statistics skipped: {conditionMessage(e)}")); NULL },
    expr  = {

      if (is_bdiv) {

        if (!inherits(params$biom, 'rbiom') || is.null(params$bdiv))
          cli_abort("Beta diversity statistics need the rbiom object (`biom`) and `bdiv` in the plot parameters.")

        bdiv_stats(
          biom         = params$biom,
          regr         = params$.xcol,
          stat.by      = params$stat.by,
          bdiv         = params$bdiv,
          tree         = params$tree,
          split.by     = setdiff(params$facet.by, c('.bdiv', '.adiv', '.rank', '.taxa')),
          strata       = params$strata,
          test         = "mantel",
          seed         = params$seed %||% 0,
          permutations = params$permutations %||% 999,
          p.adj        = params$p.adj,
          alpha        = params$alpha %||% 0.5,
          cpus         = params$cpus %||% n_cpus() )

      } else {

        stats_table(
          df       = params$.ggdata,
          regr     = params$.xcol,
          resp     = params$.ycol,
          stat.by  = params$stat.by,
          split.by = params$facet.by,
          strata   = params$strata,
          test     = params$test,
          fit      = params$fit,
          at       = params$at,
          level    = params$level,
          alt      = params$alt,
          mu       = params$mu,
          p.adj    = params$p.adj )
      }
    })

  if (is.null(stats))
    return (params)

  # Match the plot's factor levels for faceting.
  if (hasName(stats, '.bdiv') && is.factor(params$.ggdata[['.bdiv']]))
    stats[['.bdiv']] <- factor(as.character(stats[['.bdiv']]), levels = levels(params$.ggdata[['.bdiv']]))

  params$.plot_attrs[['stats']] <- stats


  #________________________________________________________
  # Use p.top to retain only most significant taxa.
  #________________________________________________________
  apply_p_top(params)
  if (is_empty_df(params$.ggdata))
    return (invisible(params))



  #________________________________________________________
  # Show p-values on the plot. One per facet (and stat.by group for
  # Mantel tests), at `at` when the test was evaluated there.
  #________________________________________________________
  xcol     <- params$.xcol
  ycol     <- params$.ycol
  stat.by  <- params$stat.by
  facet.by <- params$facet.by
  ggdata   <- params$.ggdata

  stats  <- params$.plot_attrs[['stats']]
  stats  <- stats[!is.na(stats[['.p.val']]), , drop = FALSE]

  if (nrow(stats) == 0)
    return (invisible(params))
  at_col <- !is_bdiv && hasName(stats, xcol)
  by     <- intersect(c(facet.by, if (is_bdiv) stat.by), names(stats))

  attr(ggdata, 'stat_labels') <- split_apply(stats, by, function (df) {

    df     <- df[order(df[['.p.val']]), , drop = FALSE]
    plabel <- p_text(df$.adj.p[[1]])

    if (is_bdiv)
      plabel <- paste0("Mantel r = ", format(df$.stat[[1]], digits = 2), "; ", plabel)

    if (!at_col)
      return (tibble(
        !!xcol  := -Inf,
        !!ycol  :=  Inf,
        .label   =  plabel,
        .hjust   = -0.1,
        .vjust   =  1.2 ))

    xpos   <- df[[xcol]][[1]]
    xlabel <- paste("x =", format(xpos, digits = 3))

    if (!is.null(stat.by) && hasName(df, stat.by))
      xlabel <- paste(df[[stat.by]][[1]], "\n", xlabel)

    tibble(
      !!xcol  := c(xpos, xpos),
      !!ycol  := c(Inf, -Inf),
      .label   = c(plabel, xlabel),
      .hjust   = c(0.5,  0.5),
      .vjust   = c(1, 0) )

  }) %>% as_tibble()


  attr(ggdata, 'stat_vline') <- NULL

  if (at_col)
    attr(ggdata, 'stat_vline') <- split_apply(stats, by, function (df) {
      df <- df[order(df[['.p.val']]), , drop = FALSE]
      tibble(!!xcol := df[[xcol]][[1]])
    }) %>% as_tibble()

  params$.ggdata <- ggdata


  set_layer(
    params = params,
    layer  = 'stats_label',
    'mapping|label' = ".label",
    'mapping|hjust' = ".hjust",
    'mapping|vjust' = ".vjust",
    'size'          = annot_size(params) )

  if ("linewidth" %in% names(ggplot2::GeomLabel$default_aes)) {
    set_layer(params, 'stats_label', 'linewidth' = NA)
  } else {
    set_layer(params, 'stats_label', 'label.size' = NA)
  }


  set_layer(
    params = params,
    layer  = 'yaxis',
    'expand' = c(0.15, 0, .15, 0) )

  if (!is_empty_df(attr(params$.ggdata, 'stat_vline', exact = TRUE)))
    set_layer(
      params = params,
      layer  = 'stats_vline',
      'mapping|xintercept' = params$.xcol,
      'linetype'           = 'dotdash',
      'alpha'              = 0.2 )



  if (isTRUE(params$caption)) {

    set_layer(params, 'labs', .overwrite = TRUE, caption = local({

      curr <- params$layers$labs$caption
      curr <- if (is.null(curr)) '' else paste0(curr, "\n")

      test <- params$test
      fit  <- params$fit

      interp <- if (is_bdiv) {
        sprintf("distances correlate with differences in %s (Mantel test, %s permutations)", xcol, params$permutations %||% 999)

      } else if (is.null(stat.by)) {
        if (test == "emmeans" && fit == "gam") "a non-flat trend (smooth term)"
        else if (fit == "gam")                 sprintf("trendline slope at x %s %s", params$alt, params$mu)
        else if (fit == "log")                 sprintf("trendline slope (vs log x) %s %s", params$alt, params$mu)
        else                                   sprintf("trendline slope %s %s", params$alt, params$mu)

      } else {
        if (test == "emmeans") "trendline means differ at x"
        else if (fit == "gam") "trendline slopes differ at x"
        else                   "trendline slopes differ"
      }

      if (!is.null(params$strata))
        interp <- sprintf("%s (random intercept per %s)", interp, params$strata)

      adj  <- p_adj_label(params$p.adj)
      meth <- if (is.null(adj) || nrow(params$.plot_attrs[['stats']]) <= 1) "" else
        sprintf("\np-values adjusted with the %s method.", adj)

      return (sprintf("%sLow p-value indicates %s.%s", curr, interp, meth))

    }))
  }


  return (params)
}
