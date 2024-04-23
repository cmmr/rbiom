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
#'        See [plot types][plots] for examples of each. Default: `"pe"`
#'        
#' @param color.by,shape.by,facet.by,limit.by   Metadata columns to use for 
#'        aesthetics and partitioning. See below for details. 
#'        Default: `NULL`
#'        
#' @param ...   Parameters for layer geoms (e.g. [ggplot2::geom_point()]). 
#'        Prefixing parameter names with a layer name ensures that a particular 
#'        parameter is passed to, and only to, that layer. For instance, 
#'        \code{point.size = 2} or \code{p.size = 2} ensures only the points 
#'        have their size set to \code{2}. Points can also be controlled with 
#'        the \code{pt.} prefix.
#'           
#'        
#' @return A \code{ggplot2} plot. \cr
#'         The computed sample coordinates and ggplot command 
#'         are available as \code{$data} and \code{$code} respectively. \cr
#'         If \code{color.by} is given, then \code{$stats} and 
#'         \code{$stats$code} are set. \cr
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
#'     bdiv_ord_plot(biom, layers="pemt", color.by="Body Site", rank="g")
#'     
bdiv_ord_plot <- function (
    biom, bdiv = "Bray-Curtis", ord = "UMAP", weighted = TRUE, layers = "petm", 
    color.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
    tree = NULL, test = "adonis2", seed = 0, permutations = 999, 
    rank = -1, taxa = 4, p.top = Inf, p.adj = "fdr", unc = "singly", caption = TRUE, ...) {
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("bdiv_ord_plot(%s)", as.args(params, 2, bdiv_ord_plot))
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_ord_plot', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    
    validate_bdiv(max = Inf)
    validate_rank(max = Inf)
    
    validate_unc()
    validate_ord(max = Inf)
    validate_bool('weighted', max = Inf)
    
    validate_meta_aes('color.by', null_ok = TRUE, aes = "color")
    validate_meta_aes('shape.by', null_ok = TRUE, aes = "shape", col_type = "cat")
    validate_meta_aes('facet.by', null_ok = TRUE, max = Inf,     col_type = "cat")
    validate_meta_aes('limit.by', null_ok = TRUE, max = Inf)
    
    sync_metadata()
    
    if (biom$n_samples < 4)
      stop("At least four samples are needed for an ordination.")
  })
  
  
  
  #________________________________________________________
  # Coordinates for ordination "points"
  #________________________________________________________
  with(params, {
    
    .ggdata <- bdiv_ord_table(
      biom         = biom, 
      bdiv         = bdiv,
      ord          = ord,
      weighted     = weighted,
      md           = ".all",
      split.by     = facet.by,
      stat.by      = names(color.by),
      tree         = tree,
      test         = test,
      seed         = seed,
      permutations = permutations,
      rank         = rank,
      taxa         = taxa,
      p.adj        = p.adj,
      p.top        = p.top,
      unc          = unc )
    
    .ggdata %<>% rename_cols('.sample' = ".label")
    
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
    choices = c( 'p' = "point", 'n' = "name",    's' = "spider", 
                 'n' = "name",  'd' = "density", 't' = "taxon",
                 'a' = "arrow", 'e' = "ellipse", 'm' = "mean" ) )
  
  
  
  #________________________________________________________
  # aes() parameters
  #________________________________________________________
  specs <- list(
    "arrow"      = c('x', 'y', 'xend',  'yend'),
    "spider"     = c('x', 'y', 'xend',  'yend',  'color'),
    "ellipse"    = c('x', 'y',                   'color'),
    "point"      = c('x', 'y', 'shape', 'fill',  'color'),
    "name"       = c('x', 'y',          'label', 'color'),
    "stats_text" = c(                   'label'),
    "taxon"      = c('x', 'y', 'size',  'label'),
    "mean"       = c('x', 'y', 'size') )
  
  
  args <- .qw(x, y, xend, yend, label, size)
  if (has_layer(params, 'color')) args[['color']] <- names(params$color.by)
  if (has_layer(params, 'color')) args[['fill']]  <- names(params$color.by)
  if (has_layer(params, 'shape')) args[['shape']] <- names(params$shape.by)
  
  for (layer in intersect(env_names(params$layers), names(specs))) {
    layerArgs <- args[intersect(specs[[layer]], names(args))]
    layerArgs <- args[setdiff(names(layerArgs), names(params$layers[[layer]]))]
    
    if (length(layerArgs) == 0) next
    
    names(layerArgs) <- paste0("mapping|", names(layerArgs))
    if (layer == "name") layerArgs[['mapping|label']] <- ".sample"
    set_layer(params, layer, layerArgs)
  }
  remove(list = c("specs", "args", "layer") %>% intersect(ls()))
  
  
  
  
  #________________________________________________________
  # Non-aes parameters
  #________________________________________________________
  if (has_layer(params, 'labs'))   set_layer(params, 'labs',   x = NULL, y = NULL)
  if (has_layer(params, 'spider')) set_layer(params, 'spider', alpha = 0.4, size = 0.75)
  if (has_layer(params, 'mean'))   set_layer(params, 'mean',   alpha = 0.5, color = "darkgray")
  
  if (has_layer(params, 'arrow'))
    set_layer(
      params = params, 
      layer  = 'arrow', 
      'alpha' = 0.4,
      'color' = "darkgray", 
      'size'  = 0.75, 
      'arrow' = arrow(ends="first", length=unit(.5,"cm")))
  
  if (has_layer(params, 'taxon'))
    set_layer(
      params = params, 
      layer  = 'taxon', 
      'show.legend'        = FALSE,
      'fill'               = alpha(c("white"), 0.8),
      'box.padding'        = 1,
      'segment.curvature'  = -0.1, 
      'segment.linetype'   = 8, 
      'max.overlaps'       = 100,
      'seed'               = 0)
  
  set_layer(
    params = params, 
    layer  = 'theme', 
    'axis.text'        = element_blank(),
    'axis.ticks'       = element_blank(),
    'panel.border'     = element_rect(color = "black", fill = FALSE, size = 1),
    'panel.grid.major' = element_blank(),
    'panel.grid.minor' = element_blank(),
    'panel.background' = element_rect(fill = "white"))
  
  
  
  #________________________________________________________
  # Scale the size of biplot means and taxon labels
  #________________________________________________________
  if (has_layer(params, 'mean') || has_layer(params, 'taxon')) {
    
    if (has_layer(params, 'taxon') && has_layer(params, 'mean')) {
      
      set_layer(
        params = params, 
        layer  = 'taxon', 
        'mapping|point.size' = ".size" )
      
      set_layer(
        params = params, 
        layer  = 'continuous_scale',
        'scale_name' = "size",
        'palette'    = as.cmd(scales::area_pal(range = c(2,5))),
        'aesthetics' = c("size", "point.size"), 
        'name'       = "Taxa Abundance", 
        'labels'     = ~ paste0(. * 100, "%") )
      
    } else if (has_layer(params, 'mean')) {
      
      set_layer(
        params = params, 
        layer  = 'scale_size', 
        'range'  = c(2, 5),
        'name'   = "Taxa Abundance",
        'labels' = ~ paste0(. * 100, "%") )
      
    } else {
      
      set_layer(
        params = params, 
        layer  = 'scale_size', 
        'range' = c(2, 5) )
    }
  }
  
  
  #________________________________________________________
  # Create the plot and add each layer with its arguments.
  # Also attach the human-readable ggplot command.
  #________________________________________________________
  fig <- params %>%
    ordination_biplot() %>%
    ordination_facets() %>%
    ordination_spider() %>%
    plot_build()
  
  
  
  
  attr(fig, 'cmd') <- cmd
  set_cache_value(cache_file, fig)
  
  return(fig)
}




#________________________________________________________
# Edit ggdata and taxa_coords for ggplot compatibility
#________________________________________________________
ordination_biplot <- function (params) {
  
  if (is.null(params$rank))
    return (invisible(params))
  
  
  ggdata      <- params$.ggdata
  taxa_coords <- attr(ggdata, 'taxa_coords', exact = TRUE)
  
  if (is_null(taxa_coords))
    return (invisible(params))
  
  
  #________________________________________________________
  # Duplicate ggdata rows for each biplot rank
  #________________________________________________________
  if (length(params$rank) == 1) {
    ggdata[['.rank']] <- params$rank
    
  } else {
    ggdata <- local({
      
      df <- plyr::ldply(params$rank, function (rank) {
        data.frame(check.names = FALSE, ggdata, '.rank' = rank) })
      
      for (i in setdiff(names(attributes(ggdata)), names(attributes(df))))
        attr(df, i) <- attr(ggdata, i, exact = TRUE)
      
      return (df)
    })
  }
  
  
  #________________________________________________________
  # Restructure colnames for use in geom_segment()
  #________________________________________________________
  
  final_cols <- c(
    params$facet.by, '.facet', '.weighted', '.bdiv', '.ord', 
    '.rank', '.taxa', '.x', '.y', '.x0', '.y0', '.p.val', '.adj.p', '.size' )
  
  taxa_coords %<>% keep_cols(final_cols)
  taxa_coords %<>% rename_cols('.taxa' = ".label", '.x0' = ".xend", '.y0' = ".yend")
  
  
  
  attr(ggdata, 'taxa_coords') <- taxa_coords
  params$.ggdata              <- ggdata
  
  # To enable %>% chaining
  return (invisible(params))
}




#________________________________________________________
# Faceting logic.
#________________________________________________________
ordination_facets <- function (params) {
  
  ggdata       <- params$.ggdata
  sample_stats <- params$.plot_attrs$stats
  taxa_coords  <- attr(ggdata, 'taxa_coords',  exact = TRUE)
  ranks        <- params$rank
  
  
  #________________________________________________________
  # Facet strip text = facet.by, metrics, stats, and rank
  #________________________________________________________
  df <- ggdata
  
  
  #________________________________________________________
  # Put facet.by metadata in first.
  #________________________________________________________
  for (i in params$facet.by) {
    
    df %<>% within({
      .facet <- bool_switch(
        test  = exists(".facet", inherits = FALSE), 
        yes   = sprintf("%s - %s", .facet, as.character(get(i))), 
        no    = as.character(get(i)) )
    })
  }
  
  if (hasName(df, ".facet"))
    df %<>% within(.facet %<>% sprintf(fmt="**%s**"))
  
  
  
  #________________________________________________________
  # Next are metrics with multiple values.
  #________________________________________________________
  var_metrics <- c()
  df[['.weight']] <- ifelse(df[['.weighted']], "Weighted", "Unweighted")
  for (i in c('.weight', '.bdiv', '.ord')) {
    df[[i]] %<>% as.character()
    if (length(unique(df[[i]])) > 1)
      var_metrics %<>% c(i)
  }
  
  if (length(var_metrics) == 0) {
    
    # Plot Title = "Weighted UniFrac PCoA (Phylum)"
    metric <- paste(collapse = " ", df[1, c('.weight', '.bdiv', '.ord')])
    if (length(ranks) == 1) metric <- sprintf("%s (%s)", metric, ranks)
    set_layer(params, 'labs', title = metric)
    
    
  } else if (length(var_metrics) == 1) {
    
    # Plot Title = "Weighted UniFrac (Phylum)"
    metric <- paste(collapse = " ", df[1, setdiff(c('.weight', '.bdiv', '.ord'), var_metrics)])
    if (length(ranks) == 1) metric <- sprintf("%s (%s)", metric, ranks)
    set_layer(params, 'labs', title = metric)
    
    # Facet Title = "Female - Saliva: PCoA (Phylum)"
    df[['.facet']] <- bool_switch(
      test = hasName(df, ".facet"), 
      yes  = sprintf("%s: %s", df[['.facet']], df[[var_metrics]]), 
      no   = sprintf("%s",     df[[var_metrics]]) )
    
    
  } else {
    
    # Plot Title = "Phylum"
    if (length(ranks) == 1) set_layer(params, 'labs', title = ranks)
    
    # Facet Title = "Female - Saliva<br>Weighted UniFrac PCoA (Phylum)"
    df[['.facet']] <- bool_switch(
      test = hasName(df, ".facet"), 
      yes  = sprintf("%s<br>%s %s %s", df[['.facet']], df[['.weight']], df[['.bdiv']], df[['.ord']]), 
      no   = sprintf("**%s %s %s**",                   df[['.weight']], df[['.bdiv']], df[['.ord']]) )
  }
  remove("var_metrics")
  
  
  #________________________________________________________
  # Append rank to facet title when showing multiple ranks.
  #________________________________________________________
  if (length(ranks) > 1)
    df[['.facet']] <- bool_switch(
      test  = hasName(df, ".facet"), 
      yes   = sprintf("%s (%s)", df[['.facet']], df[['.rank']]), 
      no    = sprintf("**%s**", df[['.rank']]) )
  
  
  #________________________________________________________
  # Last, add the statistics
  #________________________________________________________
  if (!is_null(sample_stats)) {
    
    # stats are agnostic of .ord and .rank
    stats_text <- with(
      data = plyr::join(df, sample_stats, by = c('.weighted', '.bdiv', params$facet.by)), 
      expr = sprintf(
        fmt = "*p* = %s; *stat* = %s; *z* = %s",
        format(.p.val, digits=3), 
        format(.stat,  digits=3), 
        format(.z,     digits=3) ))
    
    if (length(unique(stats_text)) == 1) {
      set_layer(params, 'labs',  subtitle = stats_text[[1]])
      set_layer(params, 'theme', plot.subtitle = element_markdown(size = 11))
      
    } else {
      df[['.facet']] <- bool_switch(
        test = hasName(df, ".facet"),
        yes  = sprintf("%s<br>%s", df[['.facet']], stats_text), 
        no   = stats_text )
    }
    
    
    # Add caption below plot describing stats method.
    if (isTRUE(params$caption)) {
      
      set_layer(
        params = params, 
        layer  = 'labs', 
        'caption' = sprintf(
          fmt = "Statistics computed with %s (%i permutations).", 
          params$test,
          params$permutations ))
      
      set_layer(params, 'theme', plot.caption = element_text(size = 9, face = "italic"))
    }
    
  }
  
  
  #________________________________________________________
  # Copy .facet column from df to taxa_coords
  #________________________________________________________
  if (!is_null(taxa_coords) && hasName(df, ".facet"))
    attr(df, 'taxa_coords') <- local({
      join_cols <- c('.weighted', '.bdiv', '.ord', '.rank', params$facet.by)
      join_cols <- intersect(join_cols, intersect(colnames(df), colnames(taxa_coords)))
      
      plyr::join(
          x     = taxa_coords, 
          y     = df[,c(join_cols, '.facet')], 
          by    = join_cols, 
          match = "first" ) %>% 
        drop_cols(join_cols) %>% 
        as_rbiom_tbl()
    })
  
  
  
  #________________________________________________________
  # These columns are now obsolete thanks to .facet / labs.
  #________________________________________________________
  df %<>% drop_cols('.weighted', '.bdiv', '.ord', '.rank', params$facet.by)
  
  
  #________________________________________________________
  # Overwrite facet.by with new .facet column.
  #________________________________________________________
  if (hasName(df, ".facet")) {
    
    df[['.facet']]  <- factor(df[['.facet']], levels = unique(df[['.facet']]))
    params$facet.by <- '.facet'
    
    set_layer(
      params = params, 
      layer  = 'facet', 
      'scales' = "free" )
    
    set_layer(
      params = params, 
      layer  = 'theme',
      'strip.background' = element_blank(),
      'strip.text'       = element_markdown(hjust=0) )
    
    if (!is.null(taxa_coords <- attr(df, 'taxa_coords', exact = TRUE))) {
      taxa_coords[['.facet']] %<>% factor(levels = levels(df[['.facet']]))
      attr(df, 'taxa_coords') <- taxa_coords
      remove("taxa_coords")
    }
    
  } else {
    params$facet.by <- NULL
  }
  
  
  #________________________________________________________
  # Make sure ggdata retains all its original attributes.
  #________________________________________________________
  for (i in names(attributes(ggdata)))
    if (is_null(attr(df, i, exact = TRUE)))
      attr(df, i) <- attr(ggdata, i, exact = TRUE)
  
  
  params$.ggdata <- df
  
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
  
  attr(ggdata, 'spider') <- plyr::ddply(
    .data      = as.data.frame(ggdata),
    .variables = ply_cols(c(params$facet.by, names(params$color.by))), 
    .fun       = function (df) {
      
      data.frame(
        check.names = FALSE,
        df,
        '.xend' = mean(df[['.x']]),
        '.yend' = mean(df[['.y']])
      )
    }) %>% drop_cols(".id") %>% as_tibble()
  
  params$.ggdata <- ggdata
  
  # To enable %>% chaining
  return (invisible(params))
}


