#' Ordinate samples and taxa on a 2D plane based on beta diversity distances.
#'
#' @inherit documentation_dist_test
#' @inherit documentation_default
#'
#' @family beta_diversity
#' @family ordination
#' @family visualization
#'
#'
#' @param layers   One or more of
#'        `c("point", "spider", "ellipse", "name", "mean", "taxon", "arrow")`.
#'        The first four are sample-centric; the last three are taxa-centric.
#'        Single letter abbreviations are also accepted. For instance,
#'        `c("point", "ellipse")` is equivalent to `c("p", "e")` and `"pe"`.
#'        Default: `"petm"`
#'
#' @param stat.by   Dataset field with the groups to compare and color by.
#'        Must be categorical. Default: `NULL`
#'
#' @param test    Permutational test for the sample groups (`stat.by`).
#'        Options are `"adonis2"` (PERMANOVA; the default), `"mrpp"`,
#'        `"betadisper"`, or `"none"`. See [distmat_stats()]. The taxa
#'        biplot always uses [vegan::envfit()].
#'
#' @param strata   Dataset field with the blocking/subject variable.
#'        Permutation tests permute within strata; regression models
#'        include a random intercept per stratum. Default: `NULL`
#'
#' @param ...   Parameters for layer geoms (e.g. [ggplot2::geom_point()]).
#'        Prefixing parameter names with a layer name ensures that a particular
#'        parameter is passed to, and only to, that layer. For instance,
#'        \code{point.size = 2} or \code{p.size = 2} ensures only the points
#'        have their size set to \code{2}. Points can also be controlled with
#'        the \code{pt.} prefix. An argument that no layer accepts is an
#'        error.
#'
#' @details
#' The \bold{ellipse} layer draws 95% normal-distribution ellipses
#' ([ggplot2::stat_ellipse()] with `type = "norm"`); groups with fewer than
#' four samples get no ellipse. The \bold{taxon} layer uses the ggrepel
#' package to avoid overlapping labels when it is installed, and plain
#' labels otherwise. Facet strips (or the subtitle) show the test
#' statistic, R-squared (PERMANOVA), and p-value from `$stats`.
#'
#' @return A \code{ggplot2} plot.
#'         The computed sample coordinates and ggplot command
#'         are available as \code{$data} and \code{$code} respectively.
#'         If \code{stat.by} is given, then \code{$stats} and
#'         \code{$stats$code} are set.
#'         If \code{rank} is given, then \code{$data$taxa_coords},
#'         \code{$taxa_stats}, and \code{$taxa_stats$code} are set.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     bdiv_ord_plot(biom, layers="pemt", stat.by="Body Site", rank="g")
#'
bdiv_ord_plot <- function (
    biom, bdiv = "bray", ord = "PCoA", layers = "petm",
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE,
    tree = NULL, test = "adonis2", seed = 0, permutations = 999,
    rank = -1, taxa = 4, p.top = Inf, p.adj = "fdr", unc = "singly",
    strata = NULL, caption = TRUE, alpha = 0.5, cpus = n_cpus(), ...) {

  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)

  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('bdiv_ord_plot', cmd_params(params))


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_ord_plot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {
    params$.fn <- 'bdiv_ord_plot()'
    p <- set_cache_value(cache_file, bdiv_ord_build(params))
  }

  return (plot_prov(p, cmd, biom))
}




#________________________________________________________
# bdiv_ord_plot() without the cache / provenance wrapper.
#________________________________________________________
bdiv_ord_build <- function (params) {

  stopifnot(is_bare_environment(params))

  params$.dv <- params$.data_var %||% "data"
  params$.data_code %<>% if.null(character())


  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {

    if (biom$n_samples < 4)
      cli_abort("At least four samples are needed for an ordination.")

    validate_bdiv(multiple = TRUE)
    validate_rank(max = Inf, null_ok = TRUE)
    validate_layers('layers', choices = c(
      'p' = "point",   'n' = "name",  's' = "spider",
      'd' = "density", 't' = "taxon", 'm' = "mean",
      'a' = "arrow",   'e' = "ellipse" ))

    validate_biom_field('stat.by',  null_ok = TRUE, col_type = "cat")
    validate_biom_field('facet.by', null_ok = TRUE, col_type = "cat", max = Inf)
    validate_biom_field('strata',   null_ok = TRUE, col_type = "cat")
    validate_var_choices('test', c("adonis2", "mrpp", "betadisper", "none"))
    validate_bool('caption')

    # The biplot is only computed when a taxa layer will show it.
    if (!any(c('taxon', 'arrow', 'mean') %in% layers))
      rank <- NULL


    #________________________________________________________
    # Coordinates for ordination "points"
    #________________________________________________________
    .ggdata <- bdiv_ord_table(
      biom         = biom,
      bdiv         = bdiv,
      ord          = ord,
      split.by     = facet.by,
      stat.by      = stat.by,
      strata       = strata,
      tree         = tree,
      test         = test,
      seed         = seed,
      permutations = permutations,
      rank         = rank,
      taxa         = taxa,
      p.adj        = p.adj,
      p.top        = p.top,
      unc          = unc,
      alpha        = alpha,
      cpus         = cpus )

    .xcol  <- ".x"
    .ycol  <- ".y"
    .xmode <- "numeric"


    #________________________________________________________
    # Move stats tables/commands attrs from `p$data` to `p`.
    #________________________________________________________
    .plot_attrs <- list(
      stats      = attr(.ggdata, 'stats',      exact = TRUE),
      taxa_stats = attr(.ggdata, 'taxa_stats', exact = TRUE) )
    attr(.ggdata, 'stats')      <- NULL
    attr(.ggdata, 'taxa_stats') <- NULL

  })


  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  init_layers(
    params  = params,
    choices = c( 'p' = "point",   'n' = "name",  's' = "spider",
                 'd' = "density", 't' = "taxon", 'm' = "mean",
                 'a' = "arrow",   'e' = "ellipse" ))

  if (is.null(params$rank))
    del_layer(params, c('taxon', 'arrow', 'mean'))

  repel <- pkg_available('ggrepel')


  #________________________________________________________
  # Scale the size of the biplot means by taxa abundance. Taxon labels
  # keep a fixed size, but are repelled from the means.
  #________________________________________________________
  if (has_layer(params, 'mean')) {

    if (has_layer(params, 'taxon') && repel) {

      set_layer(
        params = params,
        layer  = 'taxon',
        'mapping|point.size' = ".size" )

      set_layer(
        params = params,
        layer  = 'continuous_scale',
        'palette'    = as.cmd(scales::area_pal(range = c(1.5, 4))),
        'aesthetics' = c("size", "point.size"),
        'name'       = "Taxa Abundance",
        'labels'     = ~ paste0(. * 100, "%") )

    } else {

      set_layer(
        params = params,
        layer  = 'scale_size',
        'range'  = c(1.5, 4),
        'name'   = "Taxa Abundance",
        'labels' = ~ paste0(. * 100, "%") )
    }
  }


  #________________________________________________________
  # Provenance-tracked ggplot2 functions.
  #________________________________________________________
  .element_blank <- P('ggplot2::element_blank')
  .arrow         <- P('grid::arrow')
  .unit          <- P('grid::unit')
  .alpha         <- P('scales::alpha')


  #________________________________________________________
  # Default aes and non-aes parameters
  #________________________________________________________
  set_layer(
    params = params,
    layer  = 'ggplot',
    'mapping|x' = '.x',
    'mapping|y' = '.y' )


  if (has_layer(params, 'spider'))
    set_layer(
      params = params,
      layer  = 'spider',
      'alpha'        = 0.4,
      'linewidth'    = 0.75,
      'mapping|xend' = ".xend",
      'mapping|yend' = ".yend" )

  if (has_layer(params, 'name'))
    set_layer(
      params = params,
      layer  = 'name',
      'mapping|label' = ".sample" )

  if (has_layer(params, 'ellipse'))
    set_layer(params, 'ellipse', type = "norm")

  if (has_layer(params, 'mean'))
    set_layer(
      params = params,
      layer  = 'mean',
      'alpha'        = 0.5,
      'color'        = "darkgray",
      'mapping|size' = ".size" )

  if (has_layer(params, 'arrow'))
    set_layer(
      params = params,
      layer  = 'arrow',
      'alpha'        = 0.4,
      'color'        = "darkgray",
      'linewidth'    = 0.75,
      'arrow'        = .arrow(ends = "first", length = .unit(.5, "cm")),
      'mapping|xend' = ".xend",
      'mapping|yend' = ".yend" )

  if (has_layer(params, 'taxon')) {

    set_layer(
      params = params,
      layer  = 'taxon',
      'show.legend'   = FALSE,
      'size'          = annot_size(params),
      'fill'          = .alpha(c("white"), 0.8),
      'mapping|label' = ".label" )

    if (repel)
      set_layer(
        params = params,
        layer  = 'taxon',
        'box.padding'        = 1,
        'segment.curvature'  = -0.1,
        'segment.colour'     = "grey40",
        'max.overlaps'       = 100,
        'seed'               = 0 )
  }


  #________________________________________________________
  # Ordination axes have no meaningful units: a box without
  # tick marks or tick labels.
  #________________________________________________________
  set_layer(params, 'theme_rbiom', frame = "box")
  set_layer(params, 'theme', axis.text = .element_blank(), axis.ticks = .element_blank())


  #________________________________________________________
  # Create the plot and add each layer with its arguments.
  # Also attach the human-readable ggplot command.
  #________________________________________________________
  fig <- params %>%
    ordination_biplot() %>%
    ordination_facets() %>%
    ordination_spider() %>%
    ordination_ellipse() %>%
    plot_build()

  return (fig)
}




#________________________________________________________
# Edit ggdata and taxa_coords for ggplot compatibility
#________________________________________________________
ordination_biplot <- function (params) {

  ggdata      <- params$.ggdata
  taxa_coords <- attr(ggdata, 'taxa_coords', exact = TRUE)
  dv          <- params$.dv

  if (is.null(params$rank) || is_null(taxa_coords)) {
    attr(ggdata, 'taxa_coords') <- NULL
    params$.ggdata <- ggdata
    return (invisible(params))
  }


  #________________________________________________________
  # Restructure colnames for use in geom_segment()
  #________________________________________________________
  .taxa <- .x0 <- .y0 <- NULL # for CRAN check only

  taxa_coords <- dplyr::rename(as_tibble(taxa_coords), .label = .taxa, .xend = .x0, .yend = .y0)

  params$.data_code <- c(
    params$.data_code,
    "# Taxa biplot arrows (see `p$taxa_stats`).",
    sprintf("taxa_coords <- dplyr::rename(%s$taxa_coords, .label = .taxa, .xend = .x0, .yend = .y0)", dv) )


  #________________________________________________________
  # Duplicate ggdata rows for each biplot rank
  #________________________________________________________
  if (length(params$rank) > 1) {

    .rank <- NULL # for CRAN check only

    ggdata <- dplyr::bind_rows(lapply(params$rank, function (rank) dplyr::mutate(ggdata, .rank = rank)))

    params$.data_code <- c(
      params$.data_code,
      "",
      "# One copy of the sample coordinates per taxonomic rank.",
      sprintf("%s <- dplyr::bind_rows(lapply(%s, function (rank) dplyr::mutate(%s, .rank = rank)))", dv, deparse_vec(params$rank), dv) )
  }


  attr(ggdata, 'taxa_coords') <- taxa_coords
  params$.ggdata              <- ggdata

  # To enable %>% chaining
  return (invisible(params))
}




#________________________________________________________
# Pretty names for the beta diversity metrics.
#________________________________________________________
bdiv_pretty <- function (ids) {
  vapply(as.character(ids), FUN.VALUE = '', function (id) {
    tryCatch(ecodive::match_metric(id, div = 'beta')$name, error = function (e) id)
  })
}



#________________________________________________________
# One line of text summarizing a sample statistics row, e.g.
# "F = 10.7, R\u00b2 = 0.63, p = 0.001". With `name = TRUE`, prefixed
# by the name of the test ("PERMANOVA: F = ...").
#________________________________________________________
ord_stats_text <- function (row, test, name = TRUE) {

  fmt <- function (x, digits = 3) format(signif(x, digits), digits = digits)

  if (is.na(row[['.p.val']])) return ("statistics unavailable")

  p   <- p_text(row[['.p.val']])
  txt <- switch(
    EXPR       = test,
    adonis2    = sprintf("F = %s, R\u00b2 = %s, %s", fmt(row[['.stat']]), fmt(row[['.r.sqr']], 2), p),
    mrpp       = sprintf("delta = %s, z = %s, %s",    fmt(row[['.stat']]), fmt(row[['.z']]),         p),
    betadisper = sprintf("F = %s, %s",                fmt(row[['.stat']]), p),
    p )

  label <- switch(EXPR = test, adonis2 = "PERMANOVA", mrpp = "MRPP", betadisper = "Dispersion", NULL)
  if (isTRUE(name) && !is.null(label)) txt <- paste0(label, ": ", txt)

  return (txt)
}




#________________________________________________________
# Faceting logic. Each combination of facet.by fields, metric,
# ordination, and rank becomes one facet, labeled with the statistics.
# The labels are a lookup table (`facets`) joined into the data; the
# generated code contains it as a literal.
#________________________________________________________
ordination_facets <- function (params) {

  ggdata       <- params$.ggdata
  sample_stats <- params$.plot_attrs$stats
  taxa_coords  <- attr(ggdata, 'taxa_coords', exact = TRUE)
  axes         <- attr(ggdata, 'axes',        exact = TRUE)
  ranks        <- params$rank
  facet.by     <- params$facet.by
  dv           <- params$.dv

  #________________________________________________________
  # One row per facet.
  #________________________________________________________
  key_cols <- c(facet.by, '.bdiv', '.ord', if (length(ranks) > 1) '.rank')

  keys <- unique(as.data.frame(ggdata[key_cols]))
  keys <- as.data.frame(lapply(keys, as.character), stringsAsFactors = FALSE, optional = TRUE)
  rownames(keys) <- NULL

  labels <- rep("", nrow(keys))


  #________________________________________________________
  # Metadata facets first.
  #________________________________________________________
  if (length(facet.by) > 0)
    labels <- do.call(paste, c(keys[facet.by], sep = " - "))


  #________________________________________________________
  # Next are metrics with multiple values.
  #________________________________________________________
  var_metrics <- c('.bdiv', '.ord')[c(
    length(unique(keys[['.bdiv']])) > 1,
    length(unique(keys[['.ord']]))  > 1 )]

  metric_text <- paste(
    if ('.bdiv' %in% var_metrics) bdiv_pretty(keys[['.bdiv']]),
    if ('.ord'  %in% var_metrics) keys[['.ord']] )

  if (length(var_metrics) > 0)
    labels <- ifelse(nzchar(labels), paste0(labels, ": ", metric_text), metric_text)

  # Plot Title = "Bray-Curtis PCoA (Genus)" for the metrics that don't vary.
  # When both vary, the title is just the rank (if any).
  title <- paste(collapse = " ", c(
    if (!'.bdiv' %in% var_metrics) bdiv_pretty(keys[['.bdiv']][[1]]),
    if (!'.ord'  %in% var_metrics) as.character(keys[['.ord']][[1]]) ))
  if (length(ranks) == 1)
    title <- if (nzchar(title)) sprintf("%s (%s)", title, ranks) else as.character(ranks)
  if (nzchar(title)) set_layer(params, 'labs', title = title)


  #________________________________________________________
  # Append rank to facet title when showing multiple ranks.
  #________________________________________________________
  if (length(ranks) > 1)
    labels <- ifelse(nzchar(labels), sprintf("%s (%s)", labels, keys[['.rank']]), keys[['.rank']])


  #________________________________________________________
  # Last, add the statistics (agnostic of .ord and .rank).
  #________________________________________________________
  if (!is_null(sample_stats) && hasName(sample_stats, '.p.val')) {

    stats_by   <- intersect(c(facet.by, '.bdiv'), names(sample_stats))
    stats_keys <- do.call(paste, c(lapply(sample_stats[stats_by], as.character), sep = "\r"))
    facet_keys <- do.call(paste, c(keys[stats_by], sep = "\r"))

    # The caption names the test, so the labels don't need to.
    stats_text <- vapply(facet_keys, FUN.VALUE = '', function (k) {
      i <- match(k, stats_keys)
      if (is.na(i)) return ("statistics unavailable")
      ord_stats_text(sample_stats[i, ], params$test, name = !isTRUE(params$caption))
    })

    if (length(unique(stats_text)) == 1) {
      set_layer(params, 'labs',  subtitle = stats_text[[1]])

    } else {
      labels <- ifelse(nzchar(labels), paste0(labels, "\n", stats_text), stats_text)
    }


    # Add caption below plot describing stats method.
    if (isTRUE(params$caption)) {

      set_layer(
        params = params,
        layer  = 'labs',
        'caption' = sprintf(
          fmt = "Statistics computed with %s (%i permutations%s).",
          switch(params$test, adonis2 = "PERMANOVA", mrpp = "MRPP", betadisper = "betadisper", params$test),
          params$permutations,
          if (is.null(params$strata)) "" else sprintf(" within %s", params$strata) ))
    }

  }


  #________________________________________________________
  # Axis titles from the ordination (e.g. "PCoA 1 (24.3%)").
  #________________________________________________________
  for (axis in c('.x', '.y')) {

    lab <- unique(as.character(axes[['.label']][axes[['.axis']] == axis]))

    if (length(lab) != 1)
      lab <- if (length(unique(axes[['.ord']])) == 1)
        paste(as.character(axes[['.ord']][[1]]), match(axis, c('.x', '.y'))) else NULL

    if (axis == '.x') set_layer(params, 'labs', x = lab) else set_layer(params, 'labs', y = lab)
  }


  #________________________________________________________
  # Join the facet labels into the data (and the taxa coordinates).
  #________________________________________________________
  if (nrow(keys) > 1) {

    keys[['.facet']] <- factor(labels, levels = unique(labels))

    join_code <- function (var, cols) {
      sprintf("%s <- dplyr::left_join(dplyr::mutate(%s, dplyr::across(dplyr::all_of(%s), as.character)), facets, by = %s)",
        var, var, deparse_vec(cols), deparse_vec(cols))
    }

    join <- function (df, cols) {
      df <- dplyr::mutate(df, dplyr::across(dplyr::all_of(cols), as.character))
      dplyr::left_join(df, keys, by = cols)
    }

    params$.data_code <- c(
      params$.data_code,
      "",
      "# Facet labels (metric, ordination, and the statistics from `p$stats`).",
      sprintf("facets <- %s", tibble_code(lapply(as.list(keys), deparse_vec, indent = 4))),
      join_code(dv, key_cols) )

    attrs  <- attributes(ggdata)
    ggdata <- join(ggdata, key_cols)
    for (i in setdiff(names(attrs), c('names', 'row.names', 'class')))
      if (is.null(attr(ggdata, i, exact = TRUE)))
        attr(ggdata, i) <- attrs[[i]]

    if (!is.null(taxa_coords)) {
      tc_cols <- intersect(key_cols, colnames(taxa_coords))
      taxa_coords <- join(taxa_coords, tc_cols)
      attr(ggdata, 'taxa_coords') <- taxa_coords
      params$.data_code <- c(params$.data_code, join_code("taxa_coords", tc_cols))
    }

    params$facet.by <- '.facet'

    set_layer(params, 'facet', 'scales' = "free")

  } else {
    params$facet.by <- NULL
  }


  if (!is.null(taxa_coords))
    params$.derived_code[['taxa_coords']] <- ""  # defined in .data_code

  params$.ggdata <- ggdata

  plot_facets(params)

  # To enable %>% chaining
  return (invisible(params))
}




#________________________________________________________
# Construct data.frame for 'spider' lines
#________________________________________________________
ordination_spider <- function (params) {

  if (!has_layer(params, 'spider'))
    return (invisible(params))

  ggdata <- params$.ggdata
  dv     <- params$.dv
  by     <- unique(c(params$facet.by, params$stat.by))

  .x <- .y <- NULL # for CRAN check only

  spider <- as_tibble(ggdata)
  if (length(by) > 0) spider <- dplyr::group_by(spider, !!!rlang::syms(by))
  spider <- spider %>%
    dplyr::mutate(.xend = mean(.x), .yend = mean(.y)) %>%
    dplyr::ungroup()

  attr(ggdata, 'spider') <- spider

  params$.derived_code[['spider']] <- paste(collapse = "\n", c(
    "# Spider lines from each sample to its group's centroid.",
    sprintf("spider <- %s %%>%%", dv),
    if (length(by) > 0) sprintf("  dplyr::group_by(%s) %%>%%", paste(coan(by), collapse = ", ")),
    "  dplyr::mutate(.xend = mean(.x), .yend = mean(.y)) %>%",
    "  dplyr::ungroup()" ))

  params$.ggdata <- ggdata

  # To enable %>% chaining
  return (invisible(params))
}




#________________________________________________________
# Ellipses need at least four samples per group.
#________________________________________________________
ordination_ellipse <- function (params) {

  if (!has_layer(params, 'ellipse'))
    return (invisible(params))

  ggdata <- params$.ggdata
  dv     <- params$.dv
  by     <- unique(c(params$facet.by, params$stat.by))

  sizes <- if (length(by) == 0) nrow(ggdata) else table(do.call(paste, c(lapply(ggdata[by], as.character), sep = "\r")))

  if (all(sizes >= 4))
    return (invisible(params))

  if (!any(sizes >= 4)) {
    del_layer(params, 'ellipse')
    return (invisible(params))
  }

  ellipse <- as_tibble(ggdata)
  if (length(by) > 0) ellipse <- dplyr::group_by(ellipse, !!!rlang::syms(by))
  ellipse <- ellipse %>%
    dplyr::filter(dplyr::n() >= 4) %>%
    dplyr::ungroup()

  attr(ggdata, 'ellipse') <- ellipse
  attr(params$layers[['ellipse']], 'src') <- 'ellipse'

  params$.derived_code[['ellipse']] <- paste(collapse = "\n", c(
    "# Ellipses are only drawn for groups with at least four samples.",
    sprintf("ellipse_data <- %s %%>%%", dv),
    if (length(by) > 0) sprintf("  dplyr::group_by(%s) %%>%%", paste(coan(by), collapse = ", ")),
    "  dplyr::filter(dplyr::n() >= 4) %>%",
    "  dplyr::ungroup()" ))

  params$.ggdata <- ggdata

  return (invisible(params))
}
