

# See also ggpmisc package
# https://docs.r4photobiology.info/ggpmisc/articles/model-based-annotations.html

#' Visualize regression with scatterplots and trendlines.
#'
#' @inherit documentation_default
#'
#' @family visualization
#'
#'
#'
#' @param x   Dataset field with the x-axis values. Equivalent to the `regr`
#'        argument in [stats_table()]. Required.
#'
#' @param y   A numeric metadata column name to use for the y-axis.
#'        Default: `attr(df, 'response')`
#'
#' @param layers   One or more of
#'        `c("trend", "confidence", "point", "name", "residual")`. Single
#'        letter abbreviations are also accepted. For instance,
#'        `c("trend", "point")` is equivalent to `c("t", "p")` and `"tp"`.
#'        Default: `"tc"`
#'
#' @param test   Method for computing p-values: `'none'`, `'emmeans'`, or
#'        `'emtrends'`. See [stats_table()]. Default: `'emmeans'`
#'
#' @param strata   Dataset field with the blocking/subject variable.
#'        Permutation tests permute within strata; regression models
#'        include a random intercept per stratum. The plotted trendlines
#'        are descriptive fits without the strata term. Default: `NULL`
#'
#' @param ...   Additional parameters to pass along to ggplot2 functions.
#'        Prefix a parameter name with a layer name to pass it to only that
#'        layer. For instance, `p.size = 2` ensures only the points have their
#'        size set to `2`. Unprefixed parameters go to every layer whose
#'        function accepts them. An argument that no layer accepts is an
#'        error.
#'
#'
#' @return A `ggplot2` plot. The computed data points, ggplot2 command,
#'         stats table, and stats table commands are available as `$data`,
#'         `$code`, `$stats`, and `$stats$code`, respectively.
#'
#'
#' @section Aesthetics:
#'
#' All built-in color palettes are colorblind-friendly. The available
#' categorical palette names are: `"okabe"`, `"carto"`, `"r4"`,
#' `"polychrome"`, `"tol"`, `"bright"`, `"light"`,
#' `"muted"`, `"vibrant"`, `"tableau"`, `"classic"`,
#' `"alphabet"`, `"tableau20"`, `"kelly"`, and `"fishy"`.
#'
#' Shapes can be given as per base R - numbers 0 through 17 for various shapes,
#' or the decimal value of an ascii character, e.g. a-z = 65:90; A-Z = 97:122 to use
#' letters instead of shapes on the plot. Character strings may used as well.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- subset(hmp50, `Body Site` %in% c('Saliva', 'Stool'))
#'     df   <- adiv_table(rarefy(biom))
#'     stats_corrplot(df, "age", stat.by = "body")
#'     stats_corrplot(
#'       df       = df,
#'       x        = "Age",
#'       stat.by  = "Body Site",
#'       facet.by = "Sex",
#'       layers   = "trend" )

stats_corrplot <- function (
    df, x, y = attr(df, 'response'), layers = "tc",
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE,
    test = "emmeans", fit = "gam", at = NULL, level = 0.95, p.adj = "fdr",
    p.top = Inf, alt = "!=", mu = 0, strata = NULL, caption = TRUE, check = FALSE, ... ) {

  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('stats_corrplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- plot_cache_file('stats_corrplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {
    params$.fn        <- 'stats_corrplot()'
    params$.data_code <- data_origin_code(df)
    p <- set_cache_value(cache_file, corrplot_build(params))
  }

  return (plot_prov(p, cmd, data = df))
}




#' Visualize alpha diversity with scatterplots and trendlines.
#'
#' @inherit stats_corrplot
#' @inherit documentation_default
#'
#' @family alpha_diversity
#' @family visualization
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     p <- adiv_corrplot(babies, "age", stat.by = "deliv", fit = "gam")
#'
#'     p
#'
#'     p$stats
#'
#'     p$code

adiv_corrplot <- function (
    biom, x, adiv = "Shannon", layers = "tc",
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE,
    test = "emmeans", fit = "gam", at = NULL, level = 0.95, p.adj = "fdr",
    transform = "none", alt = "!=", mu = 0, strata = NULL, caption = TRUE, check = FALSE, ... ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('adiv_corrplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('adiv_corrplot', params)
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

      .dots$ci.min %<>% if.null(0)
      .fn <- 'adiv_corrplot()'
    })

    p <- set_cache_value(cache_file, corrplot_build(params))
  }

  return (plot_prov(p, cmd, biom))
}




#' Visualize beta diversity with scatterplots and trendlines.
#'
#' @inherit stats_corrplot
#' @inherit documentation_default
#'
#' @family beta_diversity
#' @family visualization
#'
#' @param test   Method for computing p-values: `'none'` or `'mantel'`
#'        (a Mantel test of the distances against the differences in `x`;
#'        see [bdiv_stats()]). The trendline is fitted to the pairwise
#'        distances for display only. Default: `'mantel'`
#'
#' @param strata   Dataset field with the blocking/subject variable.
#'        Permutation tests permute within strata; regression models
#'        include a random intercept per stratum. Default: `NULL`
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'     bdiv_corrplot(biom, "Age", stat.by = "Sex", layers = "tcp")

bdiv_corrplot <- function (
    biom, x, bdiv = "bray", layers = "tc",
    weighted = NULL, tree = NULL, within = NULL, between = NULL,
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE,
    test = "mantel", fit = "gam", at = NULL, level = 0.95, p.adj = "fdr",
    transform = "none", ties = "random", seed = 0, permutations = 999,
    alt = "!=", mu = 0, strata = NULL, caption = TRUE, check = FALSE,
    alpha = 0.5, cpus = n_cpus(), ... ) {

  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)

  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('bdiv_corrplot', cmd_params(params))


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_corrplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    with(params, {

      #________________________________________________________
      # Strip '==' and '!='. Append to within/between.
      #________________________________________________________
      validate_var_cmp(c('stat.by', 'facet.by'))

      validate_bdiv(multiple = TRUE)
      validate_var_choices('test', c('none', 'mantel'))
      validate_seed()
      validate_var_range('permutations', int = TRUE, range = c(1, Inf), n = 1)


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
        ties      = ties,
        seed      = seed,
        delta     = x,
        alpha     = alpha,
        cpus      = cpus )

      y <- attr(df, 'response', exact = TRUE)


      #________________________________________________________
      # Facet by metric when there are several; caption the
      # within/between selection.
      #________________________________________________________
      if (nlevels(df$.bdiv) > 1) facet.by %<>% c('.bdiv')

      .labs <- list(x = paste("Change in", x))

      if (length(c(within, between)) > 0 && isTRUE(caption))
        .labs$caption <- paste(sep = "\n", c(
          if (length(within)  > 0) paste0("within: ",  paste(within,  collapse = ", ")),
          if (length(between) > 0) paste0("between: ", paste(between, collapse = ", ")) ))

      .dots$ci.min %<>% if.null(0)
      .fn <- 'bdiv_corrplot()'
    })

    p <- set_cache_value(cache_file, corrplot_build(params))
  }

  return (plot_prov(p, cmd, biom))
}



#' Visualize rarefaction curves with scatterplots and trendlines.
#'
#' @inherit stats_corrplot
#' @inherit documentation_default
#'
#' @family rarefaction
#' @family visualization
#'
#'
#' @param fit   How to fit the trendline.
#'        Options are `'lm'`, `'log'`, and `'gam'`.
#'        Default: `'log'`
#'
#' @param rline   Where to draw a vertical line on the plot, intended to show
#'        a particular rarefaction depth. Set to `TRUE` to show an
#'        auto-selected rarefaction depth, `FALSE` to not show a line, or
#'        an integer for a custom position. Default: `TRUE`
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- subset(hmp50, `Body Site` %in% c('Saliva', 'Stool'))
#'     rare_corrplot(biom, stat.by = "body", adiv = c("sh", "o"), facet.by = "Sex")
#'

rare_corrplot <- function (
    biom, adiv = "Shannon", layers = "tc", rline = TRUE,
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE,
    test = "none", fit = "log", at = NULL, level = 0.95, p.adj = "fdr",
    transform = "none", alt = "!=", mu = 0, caption = TRUE,
    check = FALSE, cpus = n_cpus(), ... ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('rare_corrplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('rare_corrplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {
    params$.fn <- 'rare_corrplot()'
    p <- set_cache_value(cache_file, rare_corrplot_build(params))
  }

  return (plot_prov(p, cmd, biom))
}



#________________________________________________________
# rare_corrplot() without the cache / provenance wrapper, so that
# rare_multiplot() can build it as a panel.
#________________________________________________________
rare_corrplot_build <- function (params) {

  params$.dv <- params$.data_var %||% "data"

  with(params, {

    #________________________________________________________
    # Default rarefaction depth.
    #________________________________________________________
    validate_adiv(multiple = TRUE)
    validate_biom_field('stat.by',  null_ok = TRUE, col_type = "cat")
    validate_biom_field('facet.by', null_ok = TRUE, col_type = "cat", max = Inf)

    if (isTRUE(rline))  rline <- suggest_rarefy_depth(biom$counts)
    if (isFALSE(rline)) rline <- NULL
    validate_var_range('rline', n = 1, int = TRUE, null_ok = TRUE)


    #________________________________________________________
    # Select 10 depths and compute adiv metrics at each.
    #________________________________________________________
    .depths <- floor(seq(from = 5, to = stats::fivenum(sample_sums(biom))[[4]], length.out = 10))
    .md     <- unique(c(stat.by, facet.by))

    df <- dplyr::bind_rows(lapply(.depths, function (depth) {
      adiv_table(
        biom      = rarefy(biom, depth = depth),
        adiv      = adiv,
        md        = .md,
        transform = transform,
        cpus      = cpus )
    }))

    df <- strip_prov(as_rbiom_tbl(df))   # built here, not by a table function
    attr(df, 'response')   <- '.diversity'
    attr(df, 'resp_label') <- attr(adiv_table(rarefy(biom, depth = .depths[[1]]), adiv = adiv, md = NULL, transform = transform), 'resp_label', exact = TRUE)

    .data_code <- c(
      "# Alpha diversity at ten rarefaction depths.",
      "depths <- floor(seq(from = 5, to = stats::fivenum(sample_sums(biom))[[4]], length.out = 10))",
      sprintf("%s <- dplyr::bind_rows(lapply(depths, function (depth) {", .dv),
      sprintf("  %s", as.args(
        args = list(biom = aa(NA, display = "rarefy(biom, depth = depth)"), adiv = adiv, md = .md, transform = transform),
        fun  = adiv_table ) %>% sprintf(fmt = "adiv_table(%s)")),
      "}))" )

    remove(".depths", ".md")


    #________________________________________________________
    # Adjust facets and axis titles.
    #________________________________________________________
    if (nlevels(df$.adiv) > 1) facet.by %<>% c('.adiv')

    x <- '.depth'
    y <- '.diversity'

    .labs <- list(x = "Rarefaction Depth")
    if (nlevels(df$.adiv) > 1) .labs$y <- "Alpha Diversity"

    .dots$ci.min %<>% if.null(0)
  })

  return (corrplot_build(params))
}




#' Visualize taxa abundance with scatterplots and trendlines.
#'
#' @inherit stats_corrplot
#' @inherit documentation_default
#'
#' @family taxa_abundance
#' @family visualization
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(subset(hmp50, `Body Site` %in% c('Buccal mucosa', 'Saliva')))
#'     taxa_corrplot(biom, x = "BMI", stat.by = "Body Site", taxa = 'Streptococcus')

taxa_corrplot <- function (
    biom, x, rank = -1, layers = "tc",
    taxa = 6, lineage = FALSE, unc = 'singly', other = FALSE, p.top = Inf,
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE,
    test = "emmeans", fit = "gam", at = NULL, level = 0.95, p.adj = "fdr",
    transform = "none", ties = "random", seed = 0,
    alt = "!=", mu = 0, strata = NULL, caption = TRUE, check = FALSE, ... ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('taxa_corrplot', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_corrplot', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    with(params, {

      #________________________________________________________
      # Compute taxa abundance.
      #________________________________________________________
      df <- taxa_table(
        biom      = biom,
        rank      = rank,
        taxa      = taxa,
        lineage   = lineage,
        md        = unique(c(x, stat.by, facet.by, strata)),
        unc       = unc,
        other     = other,
        transform = transform,
        ties      = ties,
        seed      = seed )

      y <- attr(df, 'response', exact = TRUE)


      #________________________________________________________
      # Facet by taxon (and rank) when there are several.
      #________________________________________________________
      .labs <- list()

      if (nlevels(df$.rank) > 1) facet.by %<>% c('.rank')

      if (nlevels(df$.taxa) > 1) { facet.by %<>% c('.taxa')
      } else                     { .labs$title <- as.character(levels(df$.taxa)) }

      if (eq(transform, 'percent') || isTRUE(biom$relative))
        .dots$yaxis.labels %<>% if.null(as.cmd(scales::percent))

      .dots$ci.min %<>% if.null(0)
      .fn <- 'taxa_corrplot()'
    })

    p <- set_cache_value(cache_file, corrplot_build(params))
  }

  return (plot_prov(p, cmd, biom))
}
