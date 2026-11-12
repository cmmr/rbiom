

#' Visualize categorical metadata effects on numeric values.
#'
#' @inherit documentation_default
#'
#' @family visualization
#'
#'
#'
#' @param x   A categorical metadata column name to use for the x-axis. Or
#'        `NULL`, which groups all samples into a single category.
#'
#' @param y   A numeric metadata column name to use for the y-axis.
#'        Default: `attr(df, 'response')`
#'
#' @param layers   One or more of
#'        `c("bar", "box" ("x"), "violin", "dot", "strip", "crossbar", "errorbar", "linerange", "pointrange")`.
#'        Single letter abbreviations are also accepted. For instance,
#'        `c("box", "dot")` is equivalent to `c("x", "d")` and `"xd"`.
#'        Default: `"x"`
#'
#' @param test   Method for computing p-values. `'auto'` runs pairwise
#'        Mann-Whitney U tests (`'wilcox'`) when `x` and `stat.by` are the
#'        same field (or `'emmeans'` when `strata` is given), and
#'        Kruskal-Wallis tests (`'kruskal'`) across the `stat.by` groups
#'        within each `x` category otherwise. Beta diversity boxplots use
#'        PERMANOVA (`'adonis2'`), `'mrpp'`, or `'betadisper'` on the
#'        samples instead. `'none'` skips the statistics. See
#'        [stats_table()] and [bdiv_stats()]. Default: `'auto'`
#'
#' @param strata   Dataset field with the blocking/subject variable.
#'        Permutation tests permute within strata; regression models
#'        include a random intercept per stratum. Default: `NULL`
#'
#' @param ci   How to calculate min/max of the \bold{crossbar},
#'        \bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers.
#'        Options are: `"ci"` (confidence interval), `"range"`,
#'        `"sd"` (standard deviation), `"se"` (standard error), and
#'        `"mad"` (median absolute deviation).
#'        The center mark of \bold{crossbar} and \bold{pointrange} represents
#'        the mean, except for `"mad"` in which case it represents the median.
#'        Intervals are drawn as computed, including below zero, unless
#'        the y-axis is transformed: bounds below zero are drawn at zero
#'        on a `"sqrt"` or `"log1p"` axis, and bounds at or below zero
#'        are drawn at the smallest positive y value on a log axis.
#'        Default: `"ci"`
#'
#' @param ...   Additional parameters to pass along to ggplot2 functions.
#'        Prefix a parameter name with a layer name to pass it to only that
#'        layer. For instance, `d.size = 2` ensures only the points on the
#'        \bold{dot} layer have their size set to `2`. Unprefixed parameters
#'        go to every layer whose function accepts them (e.g.
#'        `scales = "free"` for the facets). An argument that no layer
#'        accepts is an error.
#'
#'
#' @return A `ggplot2` plot. The computed data points, ggplot2 command,
#'         stats table, and stats table commands are available as `$data`,
#'         `$code`, `$stats`, and `$stats$code`, respectively.
#'
#'
#' @section Aesthetics:
#'
#' The default colors are colorblind-friendly for up to 12 groups: Okabe and
#' Ito's palette (`"okabe"`) for up to 8, and Paul Tol's (`"tol"`) for 9-12.
#' Of the other palettes, `"bright"`, `"vibrant"`, `"muted"`, `"light"`, and
#' `"carto"` are also designed for color vision deficiencies. The available
#' categorical palette names are: `"okabe"`, `"carto"`, `"r4"`,
#' `"polychrome"`, `"tol"`, `"bright"`, `"light"`,
#' `"muted"`, `"vibrant"`, `"tableau"`, `"classic"`,
#' `"alphabet"`, `"tableau20"`, `"kelly"`, and `"fishy"`.
#'
#' Patterns are added using the fillpattern R package, which must be
#' installed to use `patterns`. Options are `"brick"`,
#' `"chevron"`, `"fish"`, `"grid"`, `"herringbone"`, `"hexagon"`, `"octagon"`,
#' `"rain"`, `"saw"`, `"shingle"`, `"rshingle"`, `"stripe"`, and `"wave"`,
#' optionally abbreviated and/or suffixed with modifiers. For example,
#' `"hex10_sm"` for the hexagon pattern rotated 10 degrees and shrunk by 2x.
#' See `fillpattern::fill_pattern()` for complete documentation of options.
#'
#' Shapes can be given as per base R - numbers 0 through 17 for various shapes,
#' or the decimal value of an ascii character, e.g. a-z = 65:90; A-Z = 97:122 to use
#' letters instead of shapes on the plot. Character strings may used as well.
#'
#' The \bold{dot} and \bold{strip} layers use the ggbeeswarm package when it
#' is installed, and jittered points otherwise.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     df <- adiv_table(rarefy(hmp50))
#'     stats_boxplot(df, x = "Body Site")
#'     stats_boxplot(df, x = "Sex", stat.by = "Body Site", layers = "be")

stats_boxplot <- function (
    df, x = NULL, y = attr(df, 'response'), layers = 'x',
    stat.by = x, facet.by = NULL, colors = TRUE, shapes = TRUE, patterns = FALSE,
    test = 'auto', flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr',
    p.top = Inf, outliers = NULL, xlab.angle = 'auto', p.label = 0.05,
    strata = NULL, caption = TRUE, ... ) {

  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('stats_boxplot', params)

  biom <- if (inherits(df, 'rbiom')) df


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- plot_cache_file('stats_boxplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {
    params$.fn        <- 'stats_boxplot()'
    params$.data_code <- data_origin_code(df)
    p <- set_cache_value(cache_file, boxplot_build(params))
  }

  return (plot_prov(p, cmd, biom = biom, data = df))
}



#' Visualize alpha diversity with boxplots.
#'
#' @inherit stats_boxplot
#' @inherit documentation_default
#'
#' @family alpha_diversity
#' @family visualization
#'
#' @export
#'
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     adiv_boxplot(biom, x="Body Site", stat.by="Body Site")
#'
#'     adiv_boxplot(biom, x="Sex", stat.by="Body Site", adiv=c("otu", "shan"), layers = "bld")
#'
#'     adiv_boxplot(biom, x="body", stat.by="sex", adiv=c('sha', 'sim'), flip=TRUE, layers="p")
#'
#'
#'     # Each plot object includes additional information.
#'     fig <- adiv_boxplot(biom, x="Body Site")
#'
#'     ## Computed Data Points -------------------
#'     fig$data
#'
#'     ## Statistics Table -----------------------
#'     fig$stats
#'
#'     ## ggplot2 Command ------------------------
#'     fig$code
#'

adiv_boxplot <- function (
    biom, x = NULL, adiv = "Shannon", layers = 'x',
    stat.by = x, facet.by = NULL, colors = TRUE, shapes = TRUE, patterns = FALSE,
    test = 'auto', flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr',
    outliers = NULL, xlab.angle = 'auto', p.label = 0.05,
    transform = "none", strata = NULL, caption = TRUE, ... ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('adiv_boxplot', params)
  schema_check('adiv_boxplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('adiv_boxplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    with(params, {

      #________________________________________________________
      # Compute alpha diversity.
      #________________________________________________________
      df <- adiv_table(
        biom      = biom,
        adiv      = adiv,
        md        = unique(c(x, stat.by, facet.by, strata)),
        transform = transform )

      y <- attr(df, 'response', exact = TRUE)


      #________________________________________________________
      # Facet by metric when there are several.
      #________________________________________________________
      if (nlevels(df$.adiv) > 1) facet.by %<>% c('.adiv')

      .fn <- 'adiv_boxplot()'
    })

    p <- set_cache_value(cache_file, boxplot_build(params))
  }

  return (plot_prov(p, cmd, biom))
}





#' Visualize BIOM data with boxplots.
#'
#' @inherit stats_boxplot
#' @inherit documentation_default
#'
#' @family beta_diversity
#' @family visualization
#'
#' @param test   Permutation test for comparing the groups. `'auto'` and
#'        `'adonis2'` run PERMANOVA; `'mrpp'` and `'betadisper'` are also
#'        available (see [bdiv_stats()]). The tests permute samples, not
#'        the plotted pairwise distances. `'none'` skips the statistics.
#'        Default: `'auto'`
#'
#' @param strata   Dataset field with the blocking/subject variable.
#'        Permutation tests permute within strata; regression models
#'        include a random intercept per stratum. Default: `NULL`
#'
#' @export
#'
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     bdiv_boxplot(biom, x="==Body Site", bdiv="unweighted_unifrac", stat.by="Body Site")

bdiv_boxplot <- function (
    biom, x = NULL, bdiv = "bray", layers = "x",
    weighted = NULL, tree = NULL, within = NULL, between = NULL,
    stat.by = x, facet.by = NULL, colors = TRUE, shapes = TRUE, patterns = FALSE,
    test = 'auto', flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr',
    outliers = NULL, xlab.angle = 'auto', p.label = 0.05,
    transform = "none", strata = NULL, seed = 0, permutations = 999,
    caption = TRUE, alpha = 0.5, cpus = n_cpus(), ... ) {

  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)

  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('bdiv_boxplot', cmd_params(params))
  schema_check('bdiv_boxplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_boxplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    with(params, {

      #________________________________________________________
      # Strip '==' and '!='. Append to within/between.
      #________________________________________________________
      validate_var_cmp(c('x', 'stat.by', 'facet.by'))

      validate_bdiv(multiple = TRUE)
      validate_seed()
      validate_var_range('permutations', int = TRUE, range = arg_range(ARGS$permutations), n = 1)


      #________________________________________________________
      # Compute beta diversity.
      #________________________________________________________
      df <- bdiv_table(
        biom      = biom,
        bdiv      = bdiv,
        tree      = tree,
        md        = unique(c(x, stat.by, facet.by)),
        within    = within,
        between   = between,
        transform = transform,
        delta     = NULL,
        alpha     = alpha,
        cpus      = cpus )

      y <- attr(df, 'response', exact = TRUE)


      #________________________________________________________
      # Facet by metric when there are several; caption the
      # within/between selection.
      #________________________________________________________
      if (nlevels(df$.bdiv) > 1) facet.by %<>% c('.bdiv')

      .labs <- list()

      if (length(c(within, between)) > 0 && isTRUE(caption))
        .labs$caption <- paste(sep = "\n", c(
          if (length(within)  > 0) paste0("within: ",  paste(within,  collapse = ", ")),
          if (length(between) > 0) paste0("between: ", paste(between, collapse = ", ")) ))

      .fn <- 'bdiv_boxplot()'
    })

    p <- set_cache_value(cache_file, boxplot_build(params))
  }

  return (plot_prov(p, cmd, biom))
}



#' Visualize BIOM data with boxplots.
#'
#' @inherit stats_boxplot
#' @inherit documentation_default
#'
#' @family taxa_abundance
#' @family visualization
#'
#'
#' @param x   A categorical metadata column name to use for the x-axis. Or
#'        `NULL`, which puts taxa along the x-axis. Default: `NULL`
#'
#' @param transform   Transformation to apply to the abundance values.
#'        Options are: `c("none", "rank", "log", "log1p", "sqrt", "percent")`.
#'        `"percent"` (and relative-abundance data) labels the y-axis as
#'        percentages. Default: `"none"`
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     taxa_boxplot(biom, stat.by = "Body Site", stripe = TRUE)
#'     taxa_boxplot(biom, layers = "bed", rank = c("Phylum", "Genus"), flip = TRUE)
#'     taxa_boxplot(
#'       biom    = subset(biom, `Body Site` %in% c('Saliva', 'Stool')),
#'       taxa    = 3,
#'       layers  = "ps",
#'       stat.by = "Body Site",
#'       colors  = c('Saliva' = "blue", 'Stool' = "red") )
#'

taxa_boxplot <- function (
    biom, x = NULL, rank = -1, layers = 'x',
    taxa = 6, unc = 'singly', other = FALSE, p.top = Inf,
    stat.by = x, facet.by = NULL, colors = TRUE, shapes = TRUE, patterns = FALSE,
    test = 'auto', flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr',
    outliers = NULL, xlab.angle = 'auto', p.label = 0.05,
    transform = 'none', y.transform = 'sqrt', strata = NULL, caption = TRUE, ... ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('taxa_boxplot', params)
  schema_check('taxa_boxplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_boxplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    with(params, {

      validate_var_choices('transform',   CHOICES$transform)
      validate_var_choices('y.transform', CHOICES$box.y.transform)


      #________________________________________________________
      # Compute taxa abundance.
      #________________________________________________________
      df <- taxa_table(
        biom      = biom,
        rank      = rank,
        taxa      = taxa,
        md        = unique(c(x, stat.by, facet.by, strata)),
        unc       = unc,
        other     = other,
        transform = transform )

      y <- attr(df, 'response', exact = TRUE)


      #________________________________________________________
      # Percent labels for relative abundances.
      #________________________________________________________
      .labs <- list()

      if (eq(transform, 'percent') || isTRUE(biom$relative))
        .dots$yaxis.labels %<>% if.null(as.cmd(scales::percent))


      #________________________________________________________
      # Facet by rank when there are several.
      #________________________________________________________
      if (nlevels(df$.rank) > 1) facet.by %<>% c('.rank')


      #________________________________________________________
      # Set x = '.taxa' or facet.by = '.taxa'.
      #________________________________________________________
      if (nlevels(df$.taxa) > 1) {
        if (is.null(x)) { x        <-   '.taxa'  }
        else            { facet.by %<>% c('.taxa') }
      }

      if (eq(x, '.taxa'))
        .labs$x <- if (nlevels(df$.rank) == 1) as.character(levels(df$.rank)) else "Taxa"

      .fn <- 'taxa_boxplot()'
    })

    p <- set_cache_value(cache_file, boxplot_build(params))
  }

  return (plot_prov(p, cmd, biom))
}
