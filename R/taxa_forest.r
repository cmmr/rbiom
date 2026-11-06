#' Plot differential abundance as a forest plot.
#'
#' Runs [taxa_da()] and shows each taxon's log2 fold change as a point with
#' its confidence interval. Significant taxa are colored by the direction
#' of the change; the others are gray. Each comparison (and rank and
#' `split.by` subset) gets its own panel, with the taxa ordered by fold
#' change.
#'
#' @inherit taxa_da params
#' @inherit documentation_default
#' @inherit documentation_plot_return return
#'
#' @family taxa_abundance
#' @family visualization
#'
#' @param p.top   Which taxa to show in each panel. A value >= 1 shows the
#'        `p.top` taxa with the smallest adjusted p-values; a value < 1
#'        shows the taxa with an adjusted p-value <= `p.top`; `Inf` shows
#'        every tested taxon. Default: `20`
#'
#' @param p.sig   Taxa with an adjusted p-value <= `p.sig` are significant
#'        and colored by the direction of the change. For ANCOM-BC2, a
#'        significant taxon must also pass its pseudo-count sensitivity
#'        analysis (`.robust`). Also passed to [taxa_da()], where it is
#'        ANCOM-BC2's `alpha`. Default: `0.05`
#'
#' @param colors   Colors for significant taxa that are more abundant in
#'        the first-named group of the comparison (or that increase with
#'        `regr`), for those that are less abundant (or decrease), and for
#'        taxa that aren't significant, in that order. Two colors use the
#'        default gray for the third. `TRUE` uses vermillion, blue, and
#'        gray; a palette name (e.g. `"okabe"`) uses its first two colors;
#'        `FALSE` draws everything in black. Default: `TRUE`
#'
#' @param caption   Add a caption describing the method, the intervals, and
#'        which taxa are shown. Default: `TRUE`
#'
#' @param ...   Parameters for the underlying ggplot2 functions. Prefix a
#'        parameter with a layer name to pass it to that layer only:
#'        `point.` ([ggplot2::geom_point()]), `ci.`
#'        ([ggplot2::geom_errorbar()]; e.g. `ci.width = 0` for no caps),
#'        `zero.` (the [ggplot2::geom_vline()] at zero), `facet.`, `labs.`,
#'        `xaxis.`, `yaxis.`, and `theme.`. An argument that no layer
#'        accepts is an error.
#'
#' @details
#' The statistics are those of [taxa_da()], which is called with the same
#' arguments; they are available as `$stats`, and the code that reproduces
#' the plot as `$code`. Each panel shows its own taxa (`p.top`), so the
#' y-axes are independent. To keep a taxon's position unique within its
#' panel, the plotted y values are the row number and taxon name (e.g.
#' `"12-Bacteroides"`); the axis labels drop the row number.
#'
#' @seealso [taxa_da()] for the statistics.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
#'     taxa_forest(biom, stat.by = "Body Site", rank = "Genus", p.top = 15)
#'
#'     biom <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool", "Buccal mucosa"))
#'     taxa_forest(biom, stat.by = "Body Site", rank = "Phylum", p.top = Inf)
#'

taxa_forest <- function (
    biom, stat.by = NULL, regr = NULL, rank = -1, taxa = NULL,
    method = "linda", covariates = NULL, split.by = NULL, strata = NULL,
    p.top = 20, p.sig = 0.05, colors = TRUE, caption = TRUE,
    prev.cut = 0.1, level = 0.95, p.adj = "fdr",
    lineage = FALSE, unc = "singly", seed = 0, ... ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('taxa_forest', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_forest', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    with(params, {
      validate_biom_field('stat.by',  null_ok = TRUE, col_type = "cat")
      validate_biom_field('regr',     null_ok = TRUE, col_type = "num")
      validate_biom_field('split.by', null_ok = TRUE, col_type = "cat", max = Inf)
      validate_rank(max = Inf)
      validate_var_range('p.top', range = c(0, Inf), n = 1)
      validate_var_range('p.sig', range = c(0, 1),   n = 1)
      validate_bool('caption')
      .fn <- 'taxa_forest()'
    })


    #________________________________________________________
    # The statistics. taxa_da() validates the remaining arguments.
    #________________________________________________________
    params$.stats <- with(params, taxa_da(
      biom       = biom,
      stat.by    = stat.by,
      regr       = regr,
      rank       = rank,
      taxa       = taxa,
      method     = method,
      covariates = covariates,
      split.by   = split.by,
      strata     = strata,
      prev.cut   = prev.cut,
      level      = level,
      p.adj      = p.adj,
      p.sig      = p.sig,
      lineage    = lineage,
      unc        = unc,
      seed       = seed ))

    p <- taxa_forest_build(params)

    set_cache_value(cache_file, p)
  }


  return (plot_prov(p, cmd, params$biom))
}




#________________________________________________________
# Select, order, and label the taxa; then draw them.
#________________________________________________________
taxa_forest_build <- function (params) {

  stats   <- params$.stats
  stat.by <- params$stat.by
  regr    <- params$regr
  p.top   <- params$p.top
  p.sig   <- params$p.sig
  robust  <- hasName(stats, '.robust')
  dv      <- "data"

  .adj.p <- .p.val <- .lfc <- NULL # for CRAN check only


  #________________________________________________________
  # Each panel is a combination of these columns.
  #________________________________________________________
  panels <- intersect(c('.rank', params$split.by, stat.by), names(stats))
  data   <- stats
  code   <- character()


  #________________________________________________________
  # The taxa to show in each panel.
  #________________________________________________________
  if (is.finite(p.top) && p.top >= 1) {

    n_top <- as.integer(floor(p.top))

    if (length(panels) > 0) data <- dplyr::group_by(data, !!!rlang::syms(panels))
    data <- dplyr::arrange(data, .adj.p, .p.val, .by_group = TRUE)
    data <- dplyr::slice_head(data, n = n_top)
    data <- dplyr::ungroup(data)

    code <- c(
      sprintf("# The %i taxa with the smallest adjusted p-values%s.", n_top, if (length(panels) > 0) " in each panel" else ""),
      sprintf("%s <- %s %%>%%", dv, dv),
      if (length(panels) > 0) sprintf("  dplyr::group_by(%s) %%>%%", paste(coan(panels), collapse = ", ")),
      "  dplyr::arrange(.adj.p, .p.val, .by_group = TRUE) %>%",
      sprintf("  dplyr::slice_head(n = %iL) %%>%%", n_top),
      "  dplyr::ungroup()",
      "" )

  } else if (is.finite(p.top)) {

    data <- dplyr::filter(data, .adj.p <= p.top)

    code <- c(
      sprintf("# The taxa with an adjusted p-value <= %s.", deparse_vec(p.top)),
      sprintf("%s <- dplyr::filter(%s, .adj.p <= %s)", dv, dv, deparse_vec(p.top)),
      "" )
  }


  #________________________________________________________
  # Order by fold change within each panel. The y values are
  # "<row>-<taxon>" so that each panel has its own positions.
  #________________________________________________________
  data       <- dplyr::arrange(data, !!!rlang::syms(c(panels, '.lfc')))
  data$.y    <- paste(seq_len(nrow(data)), data$.taxa, sep = "-")
  data$.y    <- factor(data$.y, levels = data$.y)

  code <- c(
    code,
    "# Taxa ordered by fold change; the row number keeps each panel's y positions apart.",
    sprintf("%s <- dplyr::arrange(%s, %s)", dv, dv, paste(c(coan(panels), ".lfc"), collapse = ", ")),
    sprintf("%s$.y <- paste(seq_len(nrow(%s)), %s$.taxa, sep = \"-\")", dv, dv, dv),
    sprintf("%s$.y <- factor(%s$.y, levels = %s$.y)", dv, dv, dv) )


  #________________________________________________________
  # Significance and direction of change.
  #________________________________________________________
  dirs <- forest_directions(params, stats, stat.by, regr)
  labs <- c(dirs, "Not significant")

  sig <- !is.na(data$.adj.p) & data$.adj.p <= p.sig
  if (robust) sig <- sig & data$.robust %in% TRUE

  data$.change <- ifelse(sig, ifelse(data$.lfc > 0, labs[[1]], labs[[2]]), labs[[3]])
  data$.change <- factor(data$.change, levels = intersect(labs, data$.change))

  code <- c(
    code,
    "",
    sprintf(
      "# Significant: adjusted p <= %s%s.", deparse_vec(p.sig),
      if (robust) " and robust to ANCOM-BC2's pseudo-count" else "" ),
    sprintf(
      "sig <- !is.na(%s$.adj.p) & %s$.adj.p <= %s%s", dv, dv, deparse_vec(p.sig),
      if (robust) sprintf(" & %s$.robust %%in%% TRUE", dv) else "" ),
    sprintf(
      "%s$.change <- ifelse(sig, ifelse(%s$.lfc > 0, %s, %s), %s)",
      dv, dv, double_quote(labs[[1]]), double_quote(labs[[2]]), double_quote(labs[[3]]) ),
    sprintf("%s$.change <- factor(%s$.change, levels = intersect(%s, %s$.change))", dv, dv, deparse_vec(labs), dv) )

  attr(data, 'cmd') <- attr(stats, 'cmd', exact = TRUE)

  params$.ggdata    <- data
  params$.data_code <- code
  params$.xcol      <- '.lfc'
  params$.ycol      <- '.y'
  params$.xmode     <- 'numeric'
  params$.ymode     <- 'factor'
  params$.plot_attrs <- list(stats = stats)


  #________________________________________________________
  # Colors for the significance/direction categories.
  #________________________________________________________
  params$colors   <- forest_colors(params$colors, labs, levels(data$.change))
  params$color.by <- '.change'
  params$stat.by  <- NULL


  #________________________________________________________
  # Facets: the panel columns that vary.
  #________________________________________________________
  facets <- Filter(function (col) length(unique(data[[col]])) > 1, panels)

  init_layers(params, do_init = c('zero', 'ci', 'point', if (length(facets) > 0) 'facet'))

  if (length(facets) > 0) {

    vars <- do.call(ggplot2::vars, lapply(facets, as.name))
    attr(vars, 'display') <- sprintf("vars(%s)", paste(coan(facets), collapse = ", "))
    set_layer(params, 'facet', facets = vars, scales = "free_y")

    n_panels <- nrow(unique(data[, facets, drop = FALSE]))
    dims     <- ggplot2::wrap_dims(
      n    = n_panels,
      nrow = params$layers[['facet']][['nrow']],
      ncol = params$layers[['facet']][['ncol']] )

    params$.plot_attrs$facet.nrow  <- dims[[1]]
    params$.plot_attrs$facet.ncol  <- dims[[2]]
    params$.plot_attrs$facet.count <- n_panels

  } else {
    params$.plot_attrs$facet.nrow  <- 1
    params$.plot_attrs$facet.ncol  <- 1
    params$.plot_attrs$facet.count <- 1
  }


  #________________________________________________________
  # Layers.
  #________________________________________________________
  set_layer(params, 'zero', xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3)

  set_layer(
    params    = params,
    layer     = 'ci',
    'mapping|xmin' = '.lower',
    'mapping|xmax' = '.upper',
    'mapping|y'    = '.y',
    width     = 0.3 )

  set_layer(
    params    = params,
    layer     = 'point',
    'mapping|x' = '.lfc',
    'mapping|y' = '.y',
    size      = 2 )

  set_layer(params, 'yaxis', labels = ~ sub("^[0-9]+-", "", .))
  set_layer(params, 'theme_rbiom', grid = "y")


  #________________________________________________________
  # Axis titles, legend title, and caption.
  #________________________________________________________
  ranks <- unique(as.character(if (hasName(stats, '.rank')) stats[['.rank']] else params$rank))

  set_layer(
    params = params,
    layer  = 'labs',
    x      = if (is.null(regr)) "Log2 fold change" else sprintf("Log2 fold change per unit of %s", regr),
    y      = if (length(ranks) == 1) ranks else NULL )

  if (!is.null(params$colors))
    set_layer(params, 'labs', color = sprintf("%s %s %s", forest_p_label(params$p.adj), LESS_EQUAL, format(p.sig)))

  if (isTRUE(params$caption))
    set_layer(params, 'labs', caption = forest_caption(params, stats, robust, generic = !is.null(stat.by) && eq(dirs[[1]], "Increased")))


  p <- plot_build(params)

  return (p)
}



#________________________________________________________
# Legend labels for significant increases and decreases. With a single
# two-group comparison, they name the groups.
#________________________________________________________
forest_directions <- function (params, stats, stat.by, regr) {

  if (!is.null(regr))
    return (c(sprintf("Increases with %s", regr), sprintf("Decreases with %s", regr)))

  cmps <- unique(as.character(stats[[stat.by]]))

  if (length(cmps) == 1) {
    lvls  <- levels(droplevels(as.factor(pull(params$biom, stat.by))))
    pairs <- utils::combn(lvls, 2, simplify = FALSE)
    hit   <- Filter(function (ab) eq(paste(ab[[2]], "-", ab[[1]]), cmps), pairs)
    if (length(hit) == 1)
      return (c(paste("Higher in", hit[[1]][[2]]), paste("Higher in", hit[[1]][[1]])))
  }

  c("Increased", "Decreased")
}



#________________________________________________________
# Named colors for the `.change` levels that are present.
#________________________________________________________
forest_colors <- function (colors, labs, present) {

  gray <- "#999999"

  if (is.null(colors) || isFALSE(colors)) return (NULL)
  if (length(present) == 0)               return (TRUE)

  if (isTRUE(colors)) {
    colors <- c("#D55E00", "#0072B2", gray)

  } else if (is_palette(colors)) {
    colors <- c(get_n_colors(2, colors), gray)

  } else if (is.character(colors) && length(colors) %in% 2:3) {
    colors <- c(unname(colors), gray)[1:3]

  } else {
    cli_abort("`colors` must be TRUE, FALSE, a palette name, or a character vector of two or three colors.")
  }

  stats::setNames(colors, labs)[present]
}



#________________________________________________________
# The caption: method, intervals, significance, and which taxa.
#________________________________________________________
forest_caption <- function (params, stats, robust, generic = FALSE) {

  method <- attr(stats, 'tbl_sum', exact = TRUE)[['Method']]
  p.top  <- params$p.top

  shown <- if (!is.finite(p.top)) { NULL
  } else if (p.top >= 1)          { sprintf("Showing the %i taxa with the smallest adjusted p-values per panel.", as.integer(floor(p.top)))
  } else                          { sprintf("Showing taxa with adjusted p %s %s.", LESS_EQUAL, format(p.top)) }

  paste(collapse = "\n", c(
    sprintf("%s; bars are %s%% confidence intervals.", method, format(100 * params$level)),
    if (isTRUE(generic)) "Comparisons are \"A - B\": increased taxa are more abundant in A.",
    if (robust) "Colored taxa are also robust to ANCOM-BC2's pseudo-count sensitivity analysis.",
    shown ))
}



#________________________________________________________
# The "less than or equal to" sign, kept out of the source as a literal
# so that the R code stays ASCII.
#________________________________________________________
LESS_EQUAL <- intToUtf8(0x2264)



#________________________________________________________
# "FDR-adjusted p", "Holm-adjusted p", etc.
#________________________________________________________
forest_p_label <- function (p.adj) {
  switch(
    EXPR = p.adj,
    fdr  = ,
    BH   = "FDR-adjusted p",
    BY   = "BY-adjusted p",
    none = "p",
    paste0(toupper(substr(p.adj, 1, 1)), substr(p.adj, 2, nchar(p.adj)), "-adjusted p") )
}
