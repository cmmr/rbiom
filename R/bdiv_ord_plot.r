#' Ordinate samples and taxa on a 2D plane based on beta diversity distances.
#' 
#' @name bdiv_ord_plot
#' 
#' @inherit bdiv_ord_table params
#' 
#' @family ordination
#' @family beta_diversity
#' @family visualization
#'           
#' @param layers   What graphical elements to use for drawing. Options are:
#'        \bold{point}, \bold{spider}, \bold{ellipse}, and \bold{name} for
#'        samples, and \bold{mean}, \bold{taxon}, and \bold{arrow} for taxa.
#'        Single letter abbreviations are also accepted. For instance,
#'        \code{c("point", "ellipse")} is equivalent to \code{c("p", "e")} and 
#'        \code{"pe"}.
#'        See the \code{vignette("ordination")} vignette for examples of each.
#'        Default: \code{"pe"}
#'                 
#' @param color.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for \code{vignette("data partitioning")}. Default: 
#'        \code{color.by=NULL, shape.by=NULL, facet.by=NULL, limit.by=NULL}
#'        
#' @param rank   What rank of taxa to display, for example \code{"Phylum"} or 
#'        \code{"Genus"}. Use \code{taxa_ranks()} to see all options for a 
#'        given BIOM object. The default, \code{NULL}, selects the lowest
#'        level.
#'        
#' @param ...   Parameters for layer geoms (e.g. [ggplot2::geom_point()]). 
#'        Prefixing parameter names with a layer name ensures that a particular 
#'        parameter is passed to, and only to, that layer. For instance, 
#'        \code{point.size = 2} or \code{p.size = 2} ensures only the points 
#'        have their size set to \code{2}. Points can also be controlled with 
#'        the \code{pt.} prefix.
#'        
#' @return A \code{ggplot} object.
#' 
#'         The computed sample coordinates are in \code{p$data)}.
#'         
#'         If \code{color.by} is given, then \code{attr(p, 'stats')} and 
#'         \code{attr(p, 'stats_cmds')} are set.
#'         
#'         If \code{rank} is given, then \code{attr(p$data, 'taxa_coords')}, 
#'         \code{attr(p, 'taxa_stats')}, and \code{attr(p, 'taxa_stats_cmds')} 
#'         are set.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     bdiv_ord_plot(biom, layers="pemt", color.by="Body Site", rank="g")
#'     
bdiv_ord_plot <- function (
    biom, bdiv="Bray-Curtis", ord="UMAP", weighted=TRUE, layers="pe", 
    color.by=NULL, shape.by=NULL, facet.by=NULL, limit.by=NULL, 
    tree=NULL, test="adonis2", seed=0, permutations=999, 
    rank=NULL, taxa=5, p.top=Inf, p.adj="fdr", unc="singly", ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("bdiv_ord_plot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  arg_str <- as.args(params, fun = bdiv_ord_plot, indent = 2)
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("fig  <- bdiv_ord_plot(%s)", arg_str) ))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    bdiv     %<>% validate_arg(biom, n = c(1,Inf), 'bdiv', 'bdiv', tree = tree)
    ord      %<>% validate_arg(NULL, n = c(1,Inf), 'ord')
    weighted %<>% validate_arg(NULL, n = c(1,Inf), 'weighted')
    rank     %<>% validate_arg(biom, n = c(1,Inf), 'rank', default = tail(c('OTU', taxa_ranks(biom)), 1))
  })
  
  
  #________________________________________________________
  # Subset biom by requested metadata and aes.
  #________________________________________________________
  params %<>% metadata_params(contraints = list(
    color.by   = list(n = c(0,1)),
    shape.by   = list(n = c(0,1),   col_type = "cat"),
    facet.by   = list(n = c(0,Inf), col_type = "cat"),
    limit.by   = list(n = c(0,Inf)) ))
  
  biom <- params[['biom']]
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (n_samples(biom) < 4)
    stop("At least four samples are needed for an ordination.")
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layer names.
  #________________________________________________________
  layer_names <- with(params, {
    
    layerlist <- c( 'p' = "point", 'n' = "name",    's' = "spider", 
                    'n' = "name",  'd' = "density", 't' = "taxon",
                    'a' = "arrow", 'e' = "ellipse", 'm' = "mean" )
    
    default <- if (is_null(rank)) c("p", "e") else c("t", "m")
    layer_match(layers, choices = layerlist, default = default) %>%
      c('ggplot', ., 'xaxis', 'yaxis', 'labs', 'theme', 'theme_bw')
  })
  params[['layers']] <- NULL
  
  
  #________________________________________________________
  # Disable biplot unless both rank and layers are given.
  #________________________________________________________
  if (is_null(params[['rank']]) || !any(c("taxon", "arrow", "mean") %in% layer_names)) {
    params[['rank']] <- NULL
    layer_names %<>% setdiff(c("taxon", "arrow", "mean"))
  }
  
  
  #________________________________________________________
  # Coordinates for ordination "points"
  #________________________________________________________
  ggdata <- bdiv_ord_table(
    biom         = biom, 
    bdiv         = params[['bdiv']],
    ord          = params[['ord']],
    weighted     = params[['weighted']],
    md           = TRUE,
    split.by     = params[['facet.by']],
    stat.by      = names(params[['color.by']]),
    permutations = params[['permutations']],
    rank         = params[['rank']],
    taxa         = params[['taxa']],
    p.adj        = params[['p.adj']],
    p.top        = params[['p.top']] )
  
  ggdata %<>% rename_cols('.sample' = ".label")
  
  
  
  
  #________________________________________________________
  # Initialize the `layers` object
  #________________________________________________________
  layers <- list()
  
  attr(layers, 'biom')     <- biom
  attr(layers, 'data')     <- ggdata
  attr(layers, 'params')   <- params
  attr(layers, 'function') <- bdiv_ord_plot
  attr(layers, 'xcol')     <- ".x"
  attr(layers, 'ycol')     <- ".y"
  attr(layers, 'xmode')    <- "numeric"
  attr(layers, 'ymode')    <- "numeric"
  
  if (!is_null(params[['color.by']]))   layer_names %<>% c("color")
  if (!is_null(params[['shape.by']]))   layer_names %<>% c("shape")
  if (!is_null(params[['pattern.by']])) layer_names %<>% c("pattern")
  
  initLayer(layer_names)
  # layers <- metadata_layers(layers)
  
  remove("layer_names", "ggdata")
  
  
  
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
  if (hasLayer("color")) args[['color']] <- names(params[['color.by']])
  if (hasLayer("color")) args[['fill']]  <- names(params[['color.by']])
  if (hasLayer("shape")) args[['shape']] <- names(params[['shape.by']])
  
  for (layer in intersect(names(layers), names(specs))) {
    layerArgs <- args[intersect(specs[[layer]], names(args))]
    layerArgs <- args[setdiff(names(layerArgs), names(layers[[layer]]))]
    
    if (length(layerArgs) == 0) next
    
    names(layerArgs) <- paste0("mapping|", names(layerArgs))
    setLayer(layer=layer, layerArgs)
  }
  remove(list = c("specs", "args", "arg", "layer", "val") %>% intersect(ls()))
  
  
  
  
  #________________________________________________________
  # Non-aes parameters
  #________________________________________________________
  if (hasLayer("labs"))   setLayer(layer="labs",   x = NULL, y = NULL)
  if (hasLayer("spider")) setLayer(layer="spider", alpha = 0.4, size = 0.75)
  if (hasLayer("mean"))   setLayer(layer="mean",   alpha = 0.5, color = "darkgray")
  
  if (hasLayer("arrow"))
    setLayer(layer="arrow", 
             alpha = 0.4,
             color = "darkgray", 
             size  = 0.75, 
             arrow = arrow(ends="first", length=unit(.5,"cm")))
  
  if (hasLayer("taxon"))
    setLayer(layer="taxon", 
             show.legend        = FALSE,
             fill               = alpha(c("white"), 0.8),
             box.padding        = 1,
             segment.curvature  = -0.1, 
             segment.linetype   = 8, 
             max.overlaps       = 100,
             seed               = 0)
  
  setLayer(layer="theme", 
           axis.text        = element_blank(),
           axis.ticks       = element_blank(),
           panel.border     = element_rect(color = "black", fill = FALSE, size = 1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "white"))
  
  
  
  #________________________________________________________
  # Scale the size of biplot means and taxon labels
  #________________________________________________________
  if (hasName(layers, 'mean') || hasName(layers, 'taxon')) {
    
    if (hasLayer("taxon") && hasLayer("mean")) {
      setLayer(layer="taxon", "mapping|point.size" = ".size")
      setLayer(
        layer        = "continuous_scale",
        'scale_name' = "size",
        'palette'    = as.cmd(scales::area_pal(range = c(2,5))),
        'aesthetics' = c("size", "point.size"), 
        'name'       = "Taxa Abundance", 
        'labels'     = ~ paste0(. * 100, "%") )
      
    } else if (hasLayer("mean")) {
      setLayer(
        layer    = "scale_size",
        'range'  = c(2, 5),
        'name'   = "Taxa Abundance",
        'labels' = ~ paste0(. * 100, "%") )
      
    } else {
      setLayer(layer = "scale_size", 'range' = c(2, 5))
    }
  }
  
  
  #________________________________________________________
  # Create the plot and add each layer with its arguments.
  # Also attach the human-readable ggplot command.
  #________________________________________________________
  p <- layers %>%
    ordination_biplot() %>%
    ordination_facets() %>%
    ordination_spider() %>%
    plot_build()
  
  attr(p, 'history') <- history
  
  
  #________________________________________________________
  # Move stats tables/commands attrs from `p$data` to `p`.
  #________________________________________________________
  ggdata <- p[['data']]
  
  attr(p, 'stats')           <- attr(ggdata, 'sample_stats',      exact = TRUE)
  attr(p, 'stats_cmds')      <- attr(ggdata, 'sample_stats_cmds', exact = TRUE)
  attr(p, 'taxa_stats')      <- attr(ggdata, 'taxa_stats',        exact = TRUE)
  attr(p, 'taxa_stats_cmds') <- attr(ggdata, 'taxa_stats_cmds',   exact = TRUE)
  
  attr(ggdata, 'sample_stats')      <- NULL
  attr(ggdata, 'sample_stats_cmds') <- NULL
  attr(ggdata, 'taxa_stats')        <- NULL
  attr(ggdata, 'taxa_stats_cmds')   <- NULL
  
  p[['data']] <- ggdata
  
  
  set_cache_value(cache_file, p)
  return(p)
}




#________________________________________________________
# Edit ggdata and taxa_coords for ggplot compatibility
#________________________________________________________
ordination_biplot <- function (layers) {
  
  if (!any(c("mean", "arrow", "taxon") %in% names(layers)))
    return (layers)
  
  ggdata      <- attr(layers, 'data',        exact = TRUE)
  taxa_coords <- attr(ggdata, 'taxa_coords', exact = TRUE)
  params      <- attr(layers, 'params',      exact = TRUE)
  
  if (is_null(taxa_coords))      return (layers)
  if (is_null(params[['rank']])) return (layers)
  
  
  #________________________________________________________
  # Duplicate ggdata rows for each biplot rank
  #________________________________________________________
  if (length(params[['rank']]) == 1) {
    ggdata[['.rank']] <- params[['rank']]
    
  } else {
    ggdata <- local({
      
      df <- plyr::ldply(params[['rank']], function (rank) {
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
    params[['facet.by']], '.facet', '.weighted', '.bdiv', '.ord', 
    '.rank', '.taxa', '.x', '.y', '.x0', '.y0', '.p.val', '.adj.p', '.size' )
  
  taxa_coords %<>% keep_cols(final_cols)
  taxa_coords %<>% rename_cols('.taxa' = ".label", '.x0' = ".xend", '.y0' = ".yend")
  
  
  
  attr(ggdata, 'taxa_coords') <- taxa_coords
  attr(layers, 'data')        <- ggdata
  
  return (layers)
}




#________________________________________________________
# Faceting logic.
#________________________________________________________
ordination_facets <- function (layers) {
  
  ggdata       <- attr(layers, 'data',         exact = TRUE)
  params       <- attr(layers, 'params',       exact = TRUE)
  sample_stats <- attr(ggdata, 'sample_stats', exact = TRUE)
  taxa_coords  <- attr(ggdata, 'taxa_coords',  exact = TRUE)
  ranks        <- params[['rank']]
  
  
  #________________________________________________________
  # Facet strip text = facet.by, metrics, stats, and rank
  #________________________________________________________
  df <- ggdata
  
  
  #________________________________________________________
  # Put facet.by metadata in first.
  #________________________________________________________
  for (i in params[['facet.by']]) {
    df[[i]] %<>% as.character()
    df[['.facet']] <- bool_switch(
      test  = hasName(df, ".facet"), 
      yes   = sprintf("%s - %s", df[['.facet']], df[[i]]), 
      no    = df[[i]] )
  }
  
  if (hasName(df, ".facet"))
    df[['.facet']] %<>% sprintf(fmt="**%s**")
  
  
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
    setLayer("labs", title = metric)
    
    
  } else if (length(var_metrics) == 1) {
    
    # Plot Title = "Weighted UniFrac (Phylum)"
    metric <- paste(collapse = " ", df[1, setdiff(c('.weight', '.bdiv', '.ord'), var_metrics)])
    if (length(ranks) == 1) metric <- sprintf("%s (%s)", metric, ranks)
    setLayer("labs", title = metric)
    
    # Facet Title = "Female - Saliva: PCoA (Phylum)"
    df[['.facet']] <- bool_switch(
      test = hasName(df, ".facet"), 
      yes  = sprintf("%s: %s", df[['.facet']], df[[var_metrics]]), 
      no   = sprintf("%s",     df[[var_metrics]]) )
    
    
  } else {
    
    # Plot Title = "Phylum"
    if (length(ranks) == 1) setLayer("labs", title = ranks)
    
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
    df <- plyr::join(df, sample_stats, by = c('.weighted', '.bdiv', params[['facet.by']]))
    
    stats_text <- sprintf(
      fmt = "*p* = %s; *stat* = %s; *z* = %s",
      format(df[['.p.val']], digits=3), 
      format(df[['.stat']],  digits=3), 
      format(df[['.z']],     digits=3) )
    
    methods_text <- sprintf(
      fmt = "Statistics computed with %s (%i permutations).", 
      params[['test']],
      params[['permutations']] )
    
    if (length(unique(stats_text)) == 1) {
      setLayer("labs",  subtitle = sprintf(
        fmt = "%s<br><span style='font-size:9pt'>%s</span>",
        stats_text[[1]],
        methods_text ))
      setLayer("theme", plot.subtitle = element_markdown(size = 11, lineheight = 1.2))
      
    } else {
      df[['.facet']] <- bool_switch(
        test = hasName(df, ".facet"),
        yes  = sprintf("%s<br>%s", df[['.facet']], stats_text), 
        no   = stats_text )
      setLayer("labs", caption = methods_text)
      setLayer("theme", plot.caption = element_text(face = "italic"))
    }
    
  }
  
  
  #________________________________________________________
  # Copy .facet column from df to taxa_coords
  #________________________________________________________
  if (!is_null(taxa_coords) && hasName(df, ".facet"))
    attr(ggdata, 'taxa_coords') <- local({
      join_cols <- c('.weighted', '.bdiv', '.ord', '.rank', params[['facet.by']])
      join_cols <- intersect(join_cols, intersect(colnames(df), colnames(taxa_coords)))
      
      taxa_coords <- plyr::join(
        x     = taxa_coords, 
        y     = df[,c(join_cols, '.facet'),drop=FALSE], 
        by    = join_cols, 
        match = "first" )
      
      taxa_coords %<>% drop_cols(join_cols)
      
      return (taxa_coords)
    })
  
  
  
  #________________________________________________________
  # These columns are now obsolete thanks to .facet / labs.
  #________________________________________________________
  df %<>% drop_cols('.weighted', '.bdiv', '.ord', '.rank', params[['facet.by']])
  
  
  #________________________________________________________
  # Overwrite facet.by with new .facet column.
  #________________________________________________________
  if (hasName(df, ".facet")) {
    df[['.facet']]       <- factor(df[['.facet']], levels = unique(df[['.facet']]))
    params[['facet.by']] <- '.facet'
    setLayer("facet", scales = "free")
    setLayer("theme",
             strip.background = element_blank(),
             strip.text       = element_markdown(hjust=0) )
    
    if (!is.null(taxa_coords <- attr(ggdata, 'taxa_coords', exact = TRUE))) {
      taxa_coords[['.facet']] %<>% factor(levels = levels(df[['.facet']]))
      attr(ggdata, 'taxa_coords') <- taxa_coords
      remove("taxa_coords")
    }
    
  } else {
    params[['facet.by']] <- NULL
  }
  
  
  #________________________________________________________
  # Make sure ggdata retains all its original attributes.
  #________________________________________________________
  for (i in names(attributes(ggdata)))
    if (is_null(attr(df, i, exact = TRUE)))
      attr(df, i) <- attr(ggdata, i, exact = TRUE)
  
  
  attr(layers, 'data')   <- df
  attr(layers, 'params') <- params
  
  layers <- plot_facets(layers)
  return (layers)
}




#________________________________________________________
# Constuct data.frame for 'spider' lines
#________________________________________________________
ordination_spider <- function (layers) {
  
  if (!hasName(layers, "spider")) return (layers)
  
  ggdata <- attr(layers, 'data',   exact = TRUE)
  params <- attr(layers, 'params', exact = TRUE)
  
  attr(ggdata, 'spider') <- plyr::ddply(
    .data      = ggdata,
    .variables = ply_cols(c(params[['facet.by']], names(params[['color.by']]))), 
    .fun       = function (df) {
      data.frame(
        check.names = FALSE,
        df,
        '.xend' = df[['.x']] %>% mean(),
        '.yend' = df[['.y']] %>% mean()
      )
    }) %>% drop_cols(".id")
  
  attr(layers, 'data') <- ggdata
  return (layers)
}


