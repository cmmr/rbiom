

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
#'        `'emtrends'`. Default: `'emmeans'`
#'        
#' @param ...   Additional parameters to pass along to ggplot2 functions. 
#'        Prefix a parameter name with a layer name to pass it to only that
#'        layer. For instance, `p.size = 2` ensures only the points have their 
#'        size set to `2`.
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
    p.top = Inf, alt = "!=", mu = 0, caption = TRUE, check = FALSE, ... ) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env(..., .dots = TRUE)
  cache_file <- get_cache_file('stats_corrplot', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Validate and pre-process user's arguments.
  #________________________________________________________
  params <- list2env(params)
  with(params, {
    
    if (!inherits(df, 'data.frame'))
      cli_abort("`df` must be a data.frame, not {.type {df}}.")
    
    # Enforce unnamed vectors.
    for (.i in colnames(df)) df[[.i]] %<>% unname()
    remove(".i")
    
    validate_df_field('x',        col_type = "num")
    validate_df_field('y',        col_type = "num")
    validate_df_field('stat.by',  col_type = "cat", null_ok = TRUE)
    validate_df_field('facet.by', col_type = "cat", null_ok = TRUE, max = Inf)
    
    validate_var_choices('test',  choices = c('none', 'emmeans', 'emtrends'))
    validate_var_choices('fit',   choices = c('lm', 'log', 'gam'))
    validate_var_choices('p.adj', choices = p.adjust.methods)
    validate_var_choices('alt',   choices = c("!=", ">", "<"))
    
    validate_var_range('at',    n = 1, null_ok = TRUE)
    validate_var_range('level', n = 1, range = c(0.5, 1))
    validate_var_range('mu',    n = 1)
    
    validate_bool('caption')
    
    if (!is.null(.dots$rline)) {
      rline <- .dots$rline
      .dots$rline <- NULL
    }
    
  })
  
  
  #________________________________________________________
  # Generate the plot.
  #________________________________________________________
  p <- corrplot_build(params)
  
  
  attr(p, 'cmd') <- current_cmd('stats_corrplot')
  set_cache_value(cache_file, p)
  
  
  return (p)
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
    transform = "none", alt = "!=", mu = 0, caption = TRUE, check = FALSE, ... ) {
  
  
  p <- with(slurp_env(...), {
    
    #________________________________________________________
    # Compute alpha diversity.
    #________________________________________________________
    df <- adiv_table(
      biom  = biom, 
      adiv  = adiv, 
      md    = c(x, stat.by, facet.by), 
      transform = transform )
    
    remove("biom", "adiv", "transform")
    
    
    #________________________________________________________
    # Adjust facets and y-axis title.
    #________________________________________________________
    if (nlevels(df$.adiv) > 1) facet.by %<>% c('.adiv')
    default('labs.y', attr(df, 'resp_label'))
    
    
    #________________________________________________________
    # Generate the plot using generic function.
    #________________________________________________________
    ci.min <- 0
    do.call(stats_corrplot, as.list(environment()))
    
  })
  
  attr(p, 'cmd') <- current_cmd('adiv_corrplot')
  return (p)
}




#' Visualize beta diversity with scatterplots and trendlines.
#' 
#' @inherit stats_corrplot
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family visualization
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     bdiv_corrplot(biom, "Age", stat.by = "Sex", layers = "tcp")

bdiv_corrplot <- function (
    biom, x, bdiv = "Bray-Curtis", layers = "tc", 
    weighted = TRUE, tree = NULL, within = NULL, between = NULL, 
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE, 
    test = "emmeans", fit = "gam", at = NULL, level = 0.95, p.adj = "fdr", 
    transform = "none", ties = "random", seed = 0, 
    alt = "!=", mu = 0, caption = TRUE, check = FALSE, ... ) {
  
  
  p <- with(slurp_env(...), {
    
    #________________________________________________________
    # Strip '==' and '!='. Append to within/between.
    #________________________________________________________
    validate_var_cmp(c('stat.by', 'facet.by'))
    
    
    #________________________________________________________
    # Compute beta diversity.
    #________________________________________________________
    df <- bdiv_table(
      biom      = biom,
      bdiv      = bdiv,
      weighted  = weighted,
      tree      = tree, 
      md        = c(x, stat.by, facet.by),
      within    = within,
      between   = between, 
      transform = transform, 
      ties      = ties, 
      seed      = seed,
      delta     = x ) %>%
      within({
        .bdiv <- paste(ifelse(.weighted, 'Weighted', 'Unweighted'), .bdiv)
        .bdiv <- factor(.bdiv, levels = unique(.bdiv))
        remove(".weighted") })
    
    remove("biom", "bdiv", "weighted", "tree")
    remove("transform", "ties", "seed")
    
    
    #________________________________________________________
    # Adjust facets, axis titles, and caption.
    #________________________________________________________
    if (nlevels(df$.bdiv) > 1) facet.by %<>% c('.bdiv')
    
    default('labs.x', paste("Change in", x))
    default('labs.y', attr(df, 'resp_label'))
    
    if (length(c(within, between)) > 0 && isTRUE(caption))
      labs.caption <- paste(sep = "\n", c(
        if (exists('labs.caption', inherits = FALSE)) labs.caption, 
        glue("within: {within}"), glue("between: {between}") ))
    
    remove("within", "between")
    
    
    #________________________________________________________
    # Generate the plot using generic function.
    #________________________________________________________
    ci.min <- 0
    do.call(stats_corrplot, as.list(environment()))
    
  })
  
  attr(p, 'cmd') <- current_cmd('bdiv_corrplot')
  return (p)
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
    transform = "none", alt = "!=", mu = 0, caption = TRUE, check = FALSE, ... ) {
  
  
  p <- with(slurp_env(...), {
    
    #________________________________________________________
    # Default rarefaction depth.
    #________________________________________________________
    biom <- as_rbiom(biom)
    if (isTRUE(rline))  rline <- rare_suggest(biom$counts)
    if (isFALSE(rline)) rline <- NULL
    validate_var_range('rline', n = 1, int = TRUE, null_ok = FALSE)
    
    
    #________________________________________________________
    # Select 10 depths and compute adiv metrics at each.
    #________________________________________________________
    df <- local({
      
      upper <- stats::fivenum(sample_sums(biom))[[4]]
      rLvls <- floor(seq(from = 5, to = upper, length.out = 10))
      
      plyr::ldply(rLvls, .id = ".depth", function (rLvl) {
        adiv_table(
          biom      = rarefy(biom, depth = rLvl),
          adiv      = adiv,
          md        = c(stat.by, facet.by), 
          transform = transform )
      })
      
    })
    remove("biom", "adiv", "transform")
    
    
    #________________________________________________________
    # Adjust facets and axis titles.
    #________________________________________________________
    if (nlevels(df$.adiv) > 1) facet.by %<>% c('.adiv')
    
    default('labs.x', "Rarefaction Depth")
    if (!exists('labs.y', inherits = FALSE))
      labs.y <- ifelse(
        test = nlevels(df$.adiv) > 1, 
        yes  = "Alpha Diversity", 
        no   = switch(
          EXPR = levels(df$.adiv), 
          'OTUs'  = "Observed OTUs", 
          'Depth' = "Sequencing Depth", 
          paste(levels(df$.adiv), "Diversity") ))
    
    
    #________________________________________________________
    # Generate the plot using generic function.
    #________________________________________________________
    x      <- '.depth'
    y      <- '.diversity'
    ci.min <- 0
    do.call(stats_corrplot, as.list(environment()))
    
  })
  
  attr(p, 'cmd') <- current_cmd('rare_corrplot')
  return (p)
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
    taxa = 6, lineage = FALSE, unc = 'singly', other = FALSE, 
    stat.by = NULL, facet.by = NULL, colors = TRUE, shapes = TRUE, 
    test = "emmeans", fit = "gam", at = NULL, level = 0.95, p.adj = "fdr", 
    transform = "none", ties = "random", seed = 0, 
    alt = "!=", mu = 0, caption = TRUE, check = FALSE, ... ) {
  
  params <- list2env(slurp_env(...))
  
  p <- with(params, {
    
    #________________________________________________________
    # Compute taxa abundance.
    #________________________________________________________
    df <- taxa_table(
      biom      = biom, 
      rank      = rank, 
      taxa      = taxa, 
      lineage   = lineage, 
      md        = c(x, stat.by, facet.by),
      unc       = unc, 
      other     = other,
      transform = transform, 
      ties      = ties, 
      seed      = seed )
    
    remove("biom", "rank", "taxa", "lineage", "unc", "other")
    remove("transform", "ties", "seed")
    
    
    #________________________________________________________
    # Adjust facets and y-axis title.
    #________________________________________________________
    if (nlevels(df$.taxa) > 1) { facet.by %<>% c('.taxa')
    } else { default('labs.title', levels(df$.taxa)) }
    default('labs.y', attr(df, 'resp_label'))
    
    
    #________________________________________________________
    # Generate the plot using generic function.
    #________________________________________________________
    ci.min <- 0
    do.call(stats_corrplot, as.list(environment()))
    
  })
  
  
  attr(p, 'cmd') <- current_cmd('taxa_corrplot')
  return (p)
}
