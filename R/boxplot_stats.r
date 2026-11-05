
#________________________________________________________
# Computes p-values for categorical differences and places them on the
# plot as brackets (pairwise comparisons) or labels (group-wise tests).
#
# Inputs read from `params`:
#   .ggdata, .xcol, .ycol, stat.by, facet.by, p.label, p.adj, test,
#   level, ci, outliers, flip, caption, .free_x, .free_y, p.top,
#   strata (optional; a metadata field for a random intercept / blocks),
#   and, when .ycol == ".distance" (beta diversity boxplots):
#   biom, bdiv, tree, alpha, cpus, seed, permutations.
#
# Which test runs:
#   .ycol == ".distance": PERMANOVA on the samples (bdiv_stats(); `test`
#       may be "auto"/"adonis2", "mrpp", or "betadisper"). Pairwise rows
#       become brackets between the two groups' boxes (or a label when the
#       "A vs B" pair is itself an x category); the omnibus test goes in
#       the caption.
#   x == stat.by: pairwise tests between the x categories. `test` "auto"
#       is "wilcox" (Mann-Whitney), or "emmeans" when `strata` is given;
#       "emmeans"/"wilcox"/"kruskal" are honored ("kruskal" gives one
#       group-wise label per facet).
#   x != stat.by: Kruskal-Wallis across the stat.by groups within each x
#       category (one label per category). `strata` is not supported here.
#
# Outputs:
#   params$.plot_attrs$stats      - the statistics table (see stats_table()
#                                   / bdiv_stats()).
#   attr(.ggdata, 'stat_labels')  - tibble: facet.by columns, .x, .y
#                                   (numeric), .label; used by the
#                                   'stats_text' layer.
#   attr(.ggdata, 'stat_brackets')- tibble: facet.by columns, .x, .xend,
#                                   .y, .yend; used by the 'brackets' layer.
#________________________________________________________
boxplot_stats <- function (params) {

  stopifnot(is_bare_environment(params))


  layers   <- params$layers
  ggdata   <- params$.ggdata
  xcol     <- params$.xcol
  ycol     <- params$.ycol
  stat.by  <- params$stat.by
  facet.by <- params$facet.by
  p.label  <- params$p.label
  p.adj    <- params$p.adj
  test     <- params$test %||% "auto"
  strata   <- params$strata


  #________________________________________________________
  # No work to do.
  #________________________________________________________
  if (is_empty_df(ggdata) || is.null(stat.by) || all(p.label == 0) || eq(test, 'none')) {
    set_layer(params, 'yaxis', expand = c(0.02, 0, 0.02, 0))
    return (invisible(params))
  }

  pairwise <- eq(xcol, stat.by)
  is_bdiv  <- eq(ycol, '.distance')


  #________________________________________________________
  # Beta diversity: permutation tests on the samples.
  #________________________________________________________
  if (is_bdiv) {

    method <- "permanova"
    stats  <- boxplot_bdiv_stats(params)

    if (is.null(stats)) {
      set_layer(params, 'yaxis', expand = c(0.02, 0, 0.02, 0))
      return (invisible(params))
    }
  }


  #________________________________________________________
  # Group-wise when xcol != stat.by.
  #________________________________________________________
  else if (!pairwise) {

    if (!is.null(strata))
      cli_abort(c(
        'x' = "`strata` is not supported for group-wise comparisons (`x` != `stat.by`).",
        'i' = "Set `x = stat.by` for pairwise mixed-model comparisons." ))

    method <- "kruskal"
    stats  <- stats_table(
      df       = ggdata,
      stat.by  = stat.by,
      resp     = ycol,
      test     = "kruskal",
      split.by = c(xcol, facet.by),
      p.adj    = p.adj )
  }


  #________________________________________________________
  # Pairwise.
  #________________________________________________________
  else {

    method <- switch(
      EXPR    = test,
      auto    = if (is.null(strata)) "wilcox" else "emmeans",
      wilcox  = "wilcox",
      emmeans = "emmeans",
      kruskal = "kruskal",
      cli_abort("`test` = {.val {test}} is not available for this plot; use 'wilcox', 'emmeans', or 'kruskal'.") )

    stats <- stats_table(
      df       = ggdata,
      stat.by  = xcol,
      resp     = ycol,
      test     = method,
      split.by = facet.by,
      strata   = strata,
      p.adj    = p.adj )
  }


  params$.plot_attrs[['stats']] <- stats


  #________________________________________________________
  # Implement `p.top` for taxa_boxplot().
  #________________________________________________________
  params <- apply_p_top(params)
  ggdata <- params$.ggdata
  stats  <- params$.plot_attrs[['stats']]
  if (is_empty_df(ggdata)) return (invisible(params))


  #________________________________________________________
  # Rows to display: one per bracket/label, with the x categories
  # it connects (.cat1 == .cat2 for a single-category label) and the
  # text. Omnibus PERMANOVA rows go to the caption instead.
  #________________________________________________________
  xlvls <- if (is.null(xcol)) character() else levels(ggdata[[xcol]])
  rows  <- stats

  if (hasName(rows, '.p.val')) rows <- rows[!is.na(rows[['.p.val']]), , drop = FALSE]
  if (hasName(rows, '.adj.p')) rows <- rows[rows[['.adj.p']] <= max(p.label), , drop = FALSE]

  if (method %in% c("wilcox", "emmeans")) {

    # "A - B" labels, matched against the x categories.
    pairs <- combn(xlvls, m = 2, simplify = FALSE)
    names(pairs) <- vapply(pairs, paste, '', collapse = ' - ')

    keep <- as.character(rows[[xcol]]) %in% names(pairs)
    rows <- rows[keep, , drop = FALSE]
    rows[['.cat1']] <- vapply(pairs[as.character(rows[[xcol]])], `[[`, '', 1L)
    rows[['.cat2']] <- vapply(pairs[as.character(rows[[xcol]])], `[[`, '', 2L)

  } else if (method == "permanova") {

    # "A vs B" labels: a bracket when both are x categories, a label
    # when the pair itself is an x category (between-group boxes).
    rows <- boxplot_bdiv_rows(rows, xlvls, stat.by)

  } else if (method == "kruskal") {

    if (!pairwise && !is.null(xcol) && hasName(rows, xcol)) {
      # One test per x category: label that category.
      rows[['.cat1']] <- as.character(rows[[xcol]])
      rows[['.cat2']] <- rows[['.cat1']]
    } else {
      # One test across all the x categories: a bracket spanning them.
      rows[['.cat1']] <- if (length(xlvls) > 0) xlvls[[1]]             else NA_character_
      rows[['.cat2']] <- if (length(xlvls) > 0) xlvls[[length(xlvls)]] else NA_character_
    }
  }


  #________________________________________________________
  # No p-values to display.
  #________________________________________________________
  if (is_empty_df(rows)) {
    set_layer(params, 'yaxis', expand = c(0.02, 0, 0.02, 0))
    boxplot_stats_caption(params, method, stats)
    return (invisible(params))
  }


  #________________________________________________________
  # Format how p-values are displayed.
  #________________________________________________________
  if (length(p.label) == 1) {
    rows[['.label']] <- paste("italic(p)==", formatC(rows[['.adj.p']], format = "g", digits = 1))
  } else {
    rows[['.label']] <- ""
    for (i in p.label)
      rows[['.label']] %<>% paste0(ifelse(rows[['.adj.p']] <= i, "*", ""))
  }


  #________________________________________________________
  # The y-position needed to not overlap layer content.
  #________________________________________________________
  y_pos_fun <- local({

    y_pos_funs <- c(
        if (has_layer(params, 'violin'))     "violin" else NULL,
        if (has_layer(params, 'dot'))        "max"    else NULL,
        if (has_layer(params, 'strip'))      "max"    else NULL,
        if (has_layer(params, 'box'))        "box"    else NULL,
        if (has_layer(params, 'errorbar'))   "vline"  else NULL,
        if (has_layer(params, 'pointrange')) "vline"  else NULL,
        if (has_layer(params, 'crossbar'))   "vline"  else NULL,
        if (has_layer(params, 'linerange'))  "vline"  else NULL,
        if (has_layer(params, 'bar'))        "mean"   else NULL ) %>%
      unique() %>%
      lapply(switch,
        max    = {function (v) { max(v)                        }},
        mean   = {function (v) { mean(v)                       }},
        violin = {function (v) { max(stats::density(v)[['x']]) }},
        box    = switch(as.character(isFALSE(params$outliers)),
          'TRUE'  = {function (v) { grDevices::boxplot.stats(v)$stats[5]             }},
          'FALSE' = {function (v) { max(v)                                           }} ),
        vline  = switch(params$ci,
          'range' = {function (v) { max(v)                                           }},
          'mad'   = {function (v) { median(v) + stats::mad(v, median(v))             }},
          'sd'    = {function (v) { mean(v) + sd(v)                                  }},
          'se'    = {function (v) { mean(v) + sqrt(var(v)/length(v))                 }},
          'ci'    = {function (v) { t.test(v, conf.level = params$level)$conf.int[2] }} ))

    function (df) {

      vals <- df[[ycol]]

      y_pos_funs %>%
        lapply(function (f) tryCatch(f(vals), error = function (e) NA)) %>%
        unlist() %>%
        max(0, na.rm = TRUE) %>%
        {max(0, pretty(c(0, .), n = 5), na.rm = TRUE)}
    }

  })


  #________________________________________________________
  # Minimum positions for brackets/labels, per facet.
  #________________________________________________________
  if (!isTRUE(params$.free_y)) {

    rows[['.ypos']] <- ggdata %>%
      split_apply(unique(c(xcol, stat.by, facet.by)), function (df) data.frame(.ypos = y_pos_fun(df))) %>%
      dplyr::pull('.ypos') %>%
      max(0, na.rm = TRUE)

  } else {

    rows %<>% dplyr::left_join(
      by = facet.by,
      y  = split_apply(ggdata, facet.by, function (df) {
        split_apply(df, unique(c(xcol, stat.by)), function (d) data.frame(.ypos = y_pos_fun(d))) %>%
          dplyr::pull('.ypos') %>%
          max(0, na.rm = TRUE) %>%
          data.frame(.ypos = .) }) )
  }


  #________________________________________________________
  # Geometry, per facet.
  #________________________________________________________
  rows <- split_apply(rows, facet.by, function (z) {

    # Drop x categories that are absent from this facet.
    cats <- xlvls
    if (isTRUE(params$.free_x))
      cats <- intersect(xlvls, c(as.character(z[['.cat1']]), as.character(z[['.cat2']])))

    if (length(cats) == 0) {
      # No x categories (`x = NULL`): everything sits at x = 1.
      z[['.xmin']] <- 1
      z[['.xmax']] <- 1
    } else {
      z[['.xmin']] <- as.numeric(factor(z[['.cat1']], levels = cats))
      z[['.xmax']] <- as.numeric(factor(z[['.cat2']], levels = cats))
    }

    single <- z[['.xmin']] == z[['.xmax']]

    # Single-category labels sit just above the category; brackets are
    # stacked above everything.
    z[['.step']] <- z[['.ypos']] * .13
    z[['.tick']] <- NA_real_

    s <- which(single)
    z[['.xmin']][s] <- z[['.xmin']][s] - 0.4
    z[['.xmax']][s] <- z[['.xmax']][s] + 0.4
    z[['.ypos']][s] <- z[['.ypos']][s] * 1.08
    z[['.tick']][s] <- z[['.ypos']][s]

    b <- which(!single)
    z[['.ypos']][b] <- z[['.ypos']][b] * (if (any(single)) 1.16 else 1) + z[['.step']][b] * (seq_along(b) - 1)
    z[['.tick']][b] <- z[['.ypos']][b] - z[['.step']][b] / 4
    # Round away floating point noise so the positions read well in `$code`.
    for (i in c('.ypos', '.step', '.tick'))
      z[[i]] <- signif(z[[i]], 10)

    z[, c(facet.by, '.xmin', '.xmax', '.ypos', '.step', '.tick', '.label')[c(facet.by, '.xmin', '.xmax', '.ypos', '.step', '.tick', '.label') %in% names(z)], drop = FALSE] %>%
      dplyr::select(-dplyr::any_of(facet.by))
  })


  .xmin <- .xmax <- .step <- .ypos <- .tick <- .label <- NULL # for CRAN check only

  attr(ggdata, 'stat_labels') <- tibble(
      rows,
      .x = (rows[['.xmin']] + rows[['.xmax']]) / 2,
      .y = rows[['.ypos']] + (rows[['.step']] / 15) ) %>%
    dplyr::select(-c(.xmin, .xmax, .step, .ypos, .tick))

  attr(ggdata, 'stat_brackets') <- tibble(
      dplyr::bind_rows(rows, rows, rows),
      .x    = c(rows[['.xmin']], rows[['.xmin']], rows[['.xmax']]),
      .xend = c(rows[['.xmax']], rows[['.xmin']], rows[['.xmax']]),
      .y    = c(rows[['.ypos']], rows[['.ypos']], rows[['.ypos']]),
      .yend = c(rows[['.ypos']], rows[['.tick']], rows[['.tick']]) ) %>%
    dplyr::filter(!is.na(.data$.y) & !is.na(.data$.yend)) %>%
    dplyr::select(-c(.xmin, .xmax, .step, .ypos, .tick, .label))


  set_layer(
    params = params,
    layer  = 'brackets',
    'mapping|x'    = ".x",
    'mapping|xend' = ".xend",
    'mapping|y'    = ".y",
    'mapping|yend' = ".yend" )

  set_layer(
    params = params,
    layer  = 'stats_text',
    'size'          = 3,
    'vjust'         = ifelse(isTRUE(params$flip), 0.5, 0),
    'parse'         = length(p.label) == 1,
    'mapping|x'     = ".x",
    'mapping|y'     = ".y",
    'mapping|label' = ".label" )


  boxplot_stats_caption(params, method, stats)


  #________________________________________________________
  # Add extra space to the right for the p-values.
  #________________________________________________________
  if (isTRUE(params$flip)) { set_layer(params, 'yaxis', expand = c(0.02, 0, 0.15, 0) )
  } else                   { set_layer(params, 'yaxis', expand = c(0.02, 0, 0.08, 0) ) }


  #________________________________________________________
  # Don't label the y-axis beyond the data range.
  #________________________________________________________
  if (isTRUE(params$.free_y)) {

    set_layer(params, 'yaxis', limits = c(0, NA))

  } else {

    ymax   <- min(rows[[".ypos"]])
    breaks <- pretty(c(0, ymax), n = 10)
    breaks <- breaks[breaks < ymax]

    set_layer(
      params = params,
      layer  = 'yaxis',
      breaks       = breaks[seq_along(breaks) %% 2 == 1],
      minor_breaks = breaks[seq_along(breaks) %% 2 == 0] )

    # Hide gridlines behind the p-values.
    set_layer(
      params = params,
      layer  = 'stats_bg',
      'geom'  = "rect",
      'xmin'  = -Inf,
      'xmax'  = Inf,
      'ymin'  = signif(ymax, 3),
      'ymax'  = Inf,
      'fill'  = "white",
      'color' = NA )
  }


  params$.ggdata <- ggdata

  return (invisible(params))
}




#________________________________________________________
# PERMANOVA (etc.) for a beta diversity boxplot: bdiv_stats() on the
# plot's rbiom object. Returns NULL (with a warning) when the object
# isn't available in `params`.
#________________________________________________________
boxplot_bdiv_stats <- function (params) {

  biom     <- params$biom
  bdiv     <- params$bdiv
  xcol     <- params$.xcol
  stat.by  <- params$stat.by
  facet.by <- params$facet.by
  test     <- params$test %||% "auto"

  if (!inherits(biom, 'rbiom') || is.null(bdiv)) {
    cli_warn(c(
      '!' = "Beta diversity statistics need the rbiom object (`biom`) and `bdiv` in the plot parameters; skipping.",
      'i' = "Use `bdiv_stats()` to test these groups." ))
    return (NULL)
  }

  test <- switch(
    EXPR       = test,
    auto       = "adonis2",
    adonis2    = "adonis2",
    mrpp       = "mrpp",
    betadisper = "betadisper",
    cli_abort(c(
      'x' = "`test` = {.val {test}} is not valid for beta diversity.",
      'i' = "Pairwise distances aren't independent; use 'adonis2' (PERMANOVA), 'mrpp', or 'betadisper'." )) )

  split.by <- setdiff(c(if (!eq(xcol, stat.by)) xcol, facet.by), c('.bdiv', '.adiv', '.rank', '.taxa'))

  stats <- bdiv_stats(
    biom         = biom,
    stat.by      = stat.by,
    bdiv         = bdiv,
    tree         = params$tree,
    split.by     = if (length(split.by) > 0) split.by,
    strata       = params$strata,
    test         = test,
    seed         = params$seed %||% 0,
    permutations = params$permutations %||% 999,
    p.adj        = params$p.adj,
    alpha        = params$alpha %||% 0.5,
    cpus         = params$cpus %||% n_cpus() )

  # Match the plot's factor levels for faceting.
  if (hasName(stats, '.bdiv') && is.factor(params$.ggdata[['.bdiv']]))
    stats[['.bdiv']] <- factor(as.character(stats[['.bdiv']]), levels = levels(params$.ggdata[['.bdiv']]))

  return (stats)
}



#________________________________________________________
# Which PERMANOVA rows can be drawn, and where. Pairwise rows ("A vs B")
# connect the "A" and "B" boxes when both are x categories, or label the
# "A vs B" box when that pair is an x category. Omnibus rows are dropped
# here (they go in the caption).
#________________________________________________________
boxplot_bdiv_rows <- function (rows, xlvls, stat.by) {

  rows <- rows[rows[['.groups']] == 2, , drop = FALSE]
  if (nrow(rows) == 0) return (rows)

  # Rebuild the "A vs B" labels from the group names rather than parsing
  # them, so group names containing " vs " can't be misread.
  pairs <- if (length(xlvls) >= 2) combn(xlvls, 2, simplify = FALSE) else list()
  lookup <- stats::setNames(pairs, vapply(pairs, paste, '', collapse = ' vs '))

  labels <- as.character(rows[[stat.by]])
  cat1 <- cat2 <- rep(NA_character_, nrow(rows))

  for (i in seq_along(labels)) {
    if (labels[[i]] %in% names(lookup)) {
      cat1[[i]] <- lookup[[labels[[i]]]][[1]]
      cat2[[i]] <- lookup[[labels[[i]]]][[2]]
    } else if (labels[[i]] %in% xlvls) {
      cat1[[i]] <- cat2[[i]] <- labels[[i]]
    }
  }

  rows[['.cat1']] <- cat1
  rows[['.cat2']] <- cat2
  rows[!is.na(cat1), , drop = FALSE]
}



#________________________________________________________
# The methodology caption.
#________________________________________________________
boxplot_stats_caption <- function (params, method, stats) {

  if (!isTRUE(params$caption)) return (invisible(NULL))

  .element_text <- P('ggplot2::element_text')
  set_layer(params, 'theme', plot.caption = .element_text(size = 9, face = "italic"))

  strata <- params$strata

  text <- switch(
    EXPR      = method,
    wilcox    = "Pairwise Mann-Whitney U tests",
    emmeans   = if (is.null(strata)) "Pairwise comparisons of estimated marginal means (lm)" else
                  sprintf("Pairwise comparisons of estimated marginal means (mixed model; random intercept per %s)", strata),
    kruskal   = "Kruskal-Wallis tests",
    permanova = local({

      tst <- sub("adonis2", "PERMANOVA", params$test %||% "adonis2")
      if (tst == "auto") tst <- "PERMANOVA"

      omni <- stats[stats[['.groups']] == max(stats[['.groups']]), , drop = FALSE]

      paste0(
        if (nrow(omni) == 1) sprintf("%s p = %s; ", tst, formatC(omni[['.p.val']][[1]], format = "g", digits = 2)),
        sprintf("%s with %s permutations%s", tst, params$permutations %||% 999,
          if (is.null(strata)) "" else sprintf(" within %s", strata)) )
    }) )

  adj <- p_adj_label(params$p.adj)

  text <- if (is.null(adj) || nrow(stats) <= 1) {
    paste0(text, ".")
  } else {
    sprintf("%s; p-values adjusted with the %s method.", text, adj)
  }

  # Keep any caption set by the calling function (e.g. within/between).
  curr <- params$layers[['labs']][['caption']]
  if (!is.null(curr)) text <- paste(curr, text, sep = "
")

  set_layer(params, 'labs', .overwrite = TRUE, caption = text)

  return (invisible(NULL))
}
