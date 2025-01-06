

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
#' @param test   Method for computing p-values: `'auto'` or `'none'`. `'auto'`
#'        will choose Wilcox or Kruskal-Wallis depending on the number of 
#'        groups.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 functions. 
#'        Prefix a parameter name with a layer name to pass it to only that
#'        layer. For instance, `d.size = 2` ensures only the points on the 
#'        \bold{dot} layer have their size set to `2`.
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
#' Patterns are added using the fillpattern R package. Options are `"brick"`, 
#' `"chevron"`, `"fish"`, `"grid"`, `"herringbone"`, `"hexagon"`, `"octagon"`, 
#' `"rain"`, `"saw"`, `"shingle"`, `"rshingle"`, `"stripe"`, and `"wave"`, 
#' optionally abbreviated and/or suffixed with modifiers. For example, 
#' `"hex10_sm"` for the hexagon pattern rotated 10 degrees and shrunk by 2x.
#' See [fillpattern::fill_pattern()] for complete documentation of options.
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
#'     df <- adiv_table(rarefy(hmp50))
#'     stats_boxplot(df, x = "Body Site")
#'     stats_boxplot(df, x = "Sex", stat.by = "Body Site", layers = "be")


stats_boxplot <- function (
    df, x = NULL, y = attr(df, 'response'), layers = 'x', 
    stat.by = x, facet.by = NULL, colors = TRUE, shapes = TRUE, patterns = FALSE, 
    test = 'auto', flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr', 
    p.top = Inf, outliers = NULL, xlab.angle = 'auto', p.label = 0.05, caption = TRUE, ... ) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env(..., .dots = TRUE)
  cache_file <- get_cache_file('stats_boxplot', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Generate the plot.
  #________________________________________________________
  p <- boxplot_build(list2env(params))
  
  attr(p, 'cmd') <- current_cmd('stats_boxplot')
  set_cache_value(cache_file, p)
  
  
  return (p)
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
#'     adiv_boxplot(biom, x="body", stat.by="sex", adiv=".all", flip=TRUE, layers="p")
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
    flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr', 
    outliers = NULL, xlab.angle = 'auto', p.label = 0.05, 
    transform = "none", caption = TRUE, ... ) {
  
  
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
    do.call(stats_boxplot, as.list(environment()))
    
  })
  
  attr(p, 'cmd') <- current_cmd('adiv_boxplot')
  return (p)
  
}





#' Visualize BIOM data with boxplots.
#' 
#' @inherit stats_boxplot
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family visualization
#' 
#' @export
#' 
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     bdiv_boxplot(biom, x="==Body Site", bdiv="UniFrac", stat.by="Body Site")

bdiv_boxplot <- function (
    biom, x = NULL, bdiv = "Bray-Curtis", layers = "x", 
    weighted = TRUE, tree = NULL, within = NULL, between = NULL, 
    stat.by = x, facet.by = NULL, colors = TRUE, shapes = TRUE, patterns = FALSE, 
    flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr', 
    outliers = NULL, xlab.angle = 'auto', p.label = 0.05, 
    transform = "none", caption = TRUE, ... ) {
  
  
  p <- with(slurp_env(...), {
    
    #________________________________________________________
    # Strip '==' and '!='. Append to within/between.
    #________________________________________________________
    validate_var_cmp(c('x', 'stat.by', 'facet.by'))
    
    
    #________________________________________________________
    # Compute beta diversity.
    #________________________________________________________
    df <- bdiv_table(
      biom     = biom,
      bdiv     = bdiv,
      weighted = weighted,
      tree     = tree, 
      md       = c(x, stat.by, facet.by),
      within   = within,
      between  = between, 
      transform    = transform,
      delta    = NULL ) %>%
      within({
        .bdiv <- paste(ifelse(.weighted, 'Weighted', 'Unweighted'), .bdiv)
        .bdiv <- factor(.bdiv, levels = unique(.bdiv))
        remove(".weighted") })
    
    remove("biom", "bdiv", "weighted", "tree", "transform")
    
    
    #________________________________________________________
    # Adjust facets, axis titles, and caption.
    #________________________________________________________
    if (nlevels(df$.bdiv) > 1) facet.by %<>% c('.bdiv')
    
    if (is.null(x)) {
      default('labs.x', NULL)
    } else {
      validate_df_field('x')
      default('labs.x', paste("Change in", x))
    }
    
    default('labs.y', attr(df, 'resp_label'))
    
    if (length(c(within, between)) > 0 && isTRUE(caption))
      labs.caption <- paste(sep = "\n", c(
        if (exists('labs.caption', inherits = FALSE)) labs.caption, 
        glue("within: {within}"), glue("between: {between}") ))
    
    remove("within", "between")
    
    
    #________________________________________________________
    # Generate the plot using generic function.
    #________________________________________________________
    do.call(stats_boxplot, as.list(environment()))
    
  })
  
  attr(p, 'cmd') <- current_cmd('bdiv_boxplot')
  return (p)
  
}



# To-do: add metacoder for overlaying on a phylogenetic tree.

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
    flip = FALSE, stripe = NULL, ci = 'ci', level = 0.95, p.adj = 'fdr', 
    outliers = NULL, xlab.angle = 'auto', p.label = 0.05, 
    transform = 'none', y.transform = 'sqrt', caption = TRUE, ... ) {
  
  
  p <- with(slurp_env(...), {
    
    validate_var_choices('transform',   c('none', 'rank', 'log', 'log1p', 'sqrt', 'percent'))
    validate_var_choices('y.transform', c('none', 'log1p', 'sqrt'))
    
    
    #________________________________________________________
    # Compute taxa abundance.
    #________________________________________________________
    df <- taxa_table(
      biom      = biom, 
      rank      = rank, 
      taxa      = taxa, 
      md        = c(x, stat.by, facet.by), 
      unc       = unc, 
      other     = other,
      transform = transform )
    
    remove('biom', 'rank', 'taxa', 'unc', 'other')
    
    
    #________________________________________________________
    # Re-code transformation arguments for plot_build.
    #________________________________________________________
    if (transform == 'percent') default('y.labels', as.cmd(scales::percent))
    if (y.transform == 'none' || exists('yaxis.transform')) remove('y.transform')
    remove('transform')
    
    
    #________________________________________________________
    # Adjust facets and y-axis title.
    #________________________________________________________
    if (nlevels(df$.rank) > 1) facet.by %<>% c('.rank')
    default('labs.y', attr(df, 'resp_label'))
    
    
    #________________________________________________________
    # Set x = '.taxa' or facet.by = '.taxa'.
    #________________________________________________________
    if (nlevels(df$.taxa) > 1) {
      if (is.null(x)) { x         <-    '.taxa'  }
      else            { facet.by %<>% c('.taxa') }
    }
    
    
    #________________________________________________________
    # Generate the plot using generic function.
    #________________________________________________________
    do.call(stats_boxplot, as.list(environment()))
    
  })
  
  
  attr(p, 'cmd') <- current_cmd('taxa_boxplot')
  return (p)
  
}



