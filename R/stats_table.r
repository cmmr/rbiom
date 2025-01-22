
#' Run non-parametric statistics on a data.frame.
#' 
#' A simple interface to lower-level statistics functions, including 
#' [stats::wilcox.test()], [stats::kruskal.test()], [emmeans::emmeans()], 
#' and [emmeans::emtrends()].
#' 
#' @inherit documentation_default
#' 
#' @family stats_tables
#' 
#' @return A tibble data.frame with fields from the table below. This tibble 
#' object provides the `$code` operator to print the R code used to generate 
#' the statistics.
#' 
#' | **Field**    | **Description**                                      |
#' | ------------ | ---------------------------------------------------- |
#' | .mean        | Estimated marginal mean. See [emmeans::emmeans()].   |
#' | .mean.diff   | Difference in means.                                 |
#' | .slope       | Trendline slope. See [emmeans::emtrends()].          |
#' | .slope.diff  | Difference in slopes.                                |
#' | .h1          | Alternate hypothesis.                                |
#' | .p.val       | Probability that null hypothesis is correct.         |
#' | .adj.p       | `.p.val` after adjusting for multiple comparisons.   |
#' | .effect.size | Effect size. See [emmeans::eff_size()].              |
#' | .lower       | Confidence interval lower bound.                     |
#' | .upper       | Confidence interval upper bound.                     |
#' | .se          | Standard error.                                      |
#' | .n           | Number of samples.                                   |
#' | .df          | Degrees of freedom.                                  |
#' | .stat        | Wilcoxon or Kruskal-Wallis rank sum statistic.       |
#' | .t.ratio     | `.mean` / `.se`                                      |
#' | .r.sqr       | Percent of variation explained by the model.         |
#' | .adj.r       | `.r.sqr`, taking degrees of freedom into account.    |
#' | .aic         | Akaike Information Criterion (predictive models).    |
#' | .bic         | Bayesian Information Criterion (descriptive models). |
#' | .loglik      | Log-likelihood goodness-of-fit score.                |
#' | .fit.p       | P-value for observing this fit by chance.            |
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     df <- taxa_table(biom, rank = "Family")
#'     stats_table(df, stat.by = "Body Site")[,1:6]
#'     
#'     df <- adiv_table(biom)
#'     stats_table(df, stat.by = "Sex", split.by = "Body Site")[,1:7]

stats_table <- function (
    df, regr = NULL, resp = attr(df, 'response'), 
    stat.by = NULL, split.by = NULL, 
    test = "emmeans", fit = "gam", at = NULL, 
    level = 0.95, alt = "!=", mu = 0, p.adj = "fdr" ) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- list2env(slurp_env())
  cache_file <- get_cache_file('stats_table', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  with(params, {
    
    #________________________________________________________
    # Validate and restructure user's arguments.
    #________________________________________________________
    if (!inherits(df, 'data.frame'))
      cli_abort("`df` must be a data.frame, not {.type {df}}.")
    
    validate_df_field('regr',     col_type = "num", null_ok = TRUE)
    validate_df_field('resp',     col_type = "num")
    validate_df_field('stat.by',  col_type = "cat", null_ok = TRUE)
    validate_df_field('split.by', col_type = "cat", null_ok = TRUE, max = Inf)
    
    validate_var_choices('test',  c("wilcox", "kruskal", "emmeans", "emtrends"))
    validate_var_choices('alt',   c("!=", ">", "<"))
    validate_var_choices('fit',   c("lm", "log", "gam"))
    validate_var_choices('p.adj', p.adjust.methods)
    
    validate_var_range('at',    null_ok = TRUE)
    validate_var_range('level', n = 1, range = c(0.5, 1))
    validate_var_range('mu',    n = 1)
    
    
    #________________________________________________________
    # Warn about invalid combinations of parameters.
    #________________________________________________________
    if (as.logical(anyDuplicated(c(regr, resp, stat.by, split.by))))
      cli_abort("`regr`, `resp`, `stat.by`, `split.by` must all be unique.")
    
    if (is.null(regr) && !is.null(at))
      cli_warn('`at` is ignored when `regr` = NULL.')
    
    if (test %in% c('wilcox', 'kruskal')) {
      if (!eq(fit, "gam")) cli_warn('`fit` is ignored when `test` = "{test}".')
      if (!is.null(regr))  cli_warn('`regr` is ignored when `test` = "{test}".')
    }
    
    if (!is.null(stat.by)) {
      if (!eq(level, 0.95)) cli_warn('`level` is ignored when using `stat.by`.')
      if (!eq(alt, "!="))   cli_warn('`alt` must be "!=" when using `stat.by`.')
      if (!eq(mu, 0))       cli_warn('`mu` must be 0 when using `stat.by`.')
    }
    
  })
  
  
  
  #________________________________________________________
  # Dispatch to relevant stats_* function.
  #________________________________________________________
  stats <- switch(
    EXPR = params$test,
    'wilcox'   = do.call(stats_wilcox,   fun_params(stats_wilcox,   params)), 
    'kruskal'  = do.call(stats_kruskal,  fun_params(stats_kruskal,  params)), 
    'emmeans'  = do.call(stats_emmeans,  fun_params(stats_emmeans,  params)), 
    'emtrends' = do.call(stats_emtrends, fun_params(stats_emtrends, params)) )
  
  
  attr(stats, 'cmd') <- current_cmd('stats_table')
  set_cache_value(cache_file, stats)
  
  return (stats)
}





#' Test alpha diversity for associations with metadata.
#' 
#' A convenience wrapper for [adiv_table()] + [stats_table()].
#' 
#' @inherit documentation_default
#' @inherit stats_table return
#' 
#' @family alpha_diversity
#' @family stats_tables
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- rarefy(hmp50)
#'       
#'     adiv_stats(biom, stat.by = "Sex")[,1:6]
#'       
#'     adiv_stats(biom, stat.by = "Sex", split.by = "Body Site")[,1:6]
#'     
#'     adiv_stats(biom, stat.by = "Body Site", test = "kruskal")

adiv_stats <- function (
    biom, regr = NULL, stat.by = NULL, adiv = "Shannon", 
    split.by = NULL, transform = "none", 
    test = "emmeans", fit = "gam", at = NULL, 
    level = 0.95, alt = "!=", mu = 0, p.adj = "fdr" ) {
  
  
  #________________________________________________________
  # Compute alpha diversity values
  #________________________________________________________
  df <- adiv_table(
    biom  = biom, 
    adiv  = adiv, 
    md    = c(regr, stat.by, split.by), 
    transform = transform )
  
  if (nlevels(df$.adiv) > 1)
    split.by %<>% c('.adiv')
  
  
  #________________________________________________________
  # Run statistics
  #________________________________________________________
  stats <- stats_table(
    df       = df, 
    regr     = regr, 
    stat.by  = stat.by, 
    split.by = split.by, 
    test     = test, 
    fit      = fit, 
    at       = at, 
    level    = level, 
    alt      = alt, 
    mu       = mu, 
    p.adj    = p.adj )
  
  
  attr(stats, 'cmd')  <- current_cmd('adiv_stats')
  attr(stats, 'data') <- df
  
  return (stats)
  
}



#' Test beta diversity for associations with metadata.
#' 
#' A convenience wrapper for [bdiv_table()] + [stats_table()].
#' 
#' @inherit documentation_default
#' @inherit stats_table return
#' 
#' @family beta_diversity
#' @family stats_tables
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'       
#'     bdiv_stats(biom, stat.by = "Sex", bdiv = c("bray", "unifrac"))[,1:7]
#'     
#'     biom <- subset(biom, `Body Site` %in% c('Saliva', 'Stool', 'Buccal mucosa'))
#'     bdiv_stats(biom, stat.by = "Body Site", split.by = "==Sex")[,1:6]

bdiv_stats <- function (
    biom, regr = NULL, stat.by = NULL, bdiv = "Bray-Curtis", 
    weighted = TRUE, tree = NULL, within = NULL, between = NULL, 
    split.by = NULL, transform = "none", 
    test = "emmeans", fit = "gam", at = NULL, 
    level = 0.95, alt = "!=", mu = 0, p.adj = "fdr" ) {
  
  
  #________________________________________________________
  # Compute beta diversity values
  #________________________________________________________
  df <- bdiv_table(
    biom      = biom, 
    bdiv      = bdiv, 
    weighted  = weighted, 
    tree      = tree, 
    md        = c(regr, stat.by, split.by, within, between), 
    within    = within, 
    between   = between, 
    delta     = regr, 
    transform = transform )
  
  if (nlevels(df$.bdiv) > 1)
    split.by %<>% c('.bdiv')
  
  
  #________________________________________________________
  # Run statistics
  #________________________________________________________
  stats <- stats_table(
    df       = df, 
    regr     = regr, 
    stat.by  = if (!is.null(stat.by))  sub("^(==|!=)", "", stat.by), 
    split.by = if (!is.null(split.by)) sub("^(==|!=)", "", split.by), 
    test     = test, 
    fit      = fit, 
    at       = at, 
    level    = level, 
    alt      = alt, 
    mu       = mu, 
    p.adj    = p.adj )
  
  
  attr(stats, 'cmd') <- current_cmd('bdiv_stats')
  attr(stats, 'data') <- df
  
  return (stats)
  
}





#' Test taxa abundances for associations with metadata.
#' 
#' A convenience wrapper for [taxa_table()] + [stats_table()].
#' 
#' @inherit documentation_default
#' @inherit stats_table return
#' 
#' @family taxa_abundance
#' @family stats_tables
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     taxa_stats(biom, stat.by = "Body Site", rank = "Family")[,1:6]

taxa_stats <- function (
    biom, regr = NULL, stat.by = NULL, rank = -1, taxa = 6, 
    lineage = FALSE, unc = "singly", other = FALSE,
    split.by = NULL, transform = "none", 
    test = "emmeans", fit = "gam", at = NULL, 
    level = 0.95, alt = "!=", mu = 0, p.adj = "fdr" ) {
  
  
  #________________________________________________________
  # Compute taxa abundance values
  #________________________________________________________
  df <- taxa_table(
    biom    = biom, 
    rank    = rank, 
    taxa    = taxa, 
    lineage = lineage, 
    md      = setdiff(c(regr, stat.by, split.by), ".taxa"), 
    unc     = unc, 
    other   = other, 
    transform   = transform )
  
  if (nlevels(df$.rank) > 1) split.by %<>% c('.rank')
  if (!eq(stat.by, ".taxa")) split.by %<>% c('.taxa')
  
  
  #________________________________________________________
  # Run statistics
  #________________________________________________________
  stats <- stats_table(
    df       = df, 
    regr     = regr, 
    stat.by  = stat.by, 
    split.by = split.by, 
    test     = test, 
    fit      = fit, 
    at       = at, 
    level    = level, 
    alt      = alt, 
    mu       = mu, 
    p.adj    = p.adj )
  
  
  attr(stats, 'cmd') <- current_cmd('taxa_stats')
  attr(stats, 'data') <- df
  
  return (stats)
  
}



