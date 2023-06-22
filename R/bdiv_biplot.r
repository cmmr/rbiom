#' Display distances between samples as a 2-D scatter plot.
#' 
#' @name bdiv_biplot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param metric   Beta diversity metric to use for calculating inter-sample
#'        distances. Options are: \code{"Bray-Curtis"}, \code{"Manhattan"},
#'        \code{"Euclidean"}, \code{"Jaccard"}, or \code{"UniFrac"}.
#'        Default: \code{"Bray-Curtis"}.
#'        
#' @param ord   Ordination method to use for condensing n-dimensional data to
#'        two dimensions. Options are: \code{"UMAP"}, \code{"tSNE"},
#'        \code{"NMDS"}, or \code{"PCoA"}. Default: \code{"PCoA"}.
#'           
#' @param layers   What graphical elements to use for drawing. Options are:
#'        \bold{point}, \bold{spider}, \bold{ellipse}, and \bold{name} for
#'        samples, and \bold{mean}, \bold{taxon}, and \bold{arrow} for biplots.
#'        Single letter abbreviations are also accepted. For instance,
#'        \code{c("point", "ellipse")} is equivalent to \code{c("p", "e")} and 
#'        \code{"pe"}.
#'        See \code{vignette("biplots")} for examples of each.
#'        Default: \code{"pce"}/\code{"p"} for ordinations 
#'        with/without a \code{color.by} argument.
#'                 
#' @param color.by,shape.by,facet.by,limit.by   Metadata columns to 
#'        use for data partitioning. Default: \code{NULL}
#'        
#' @param weighted   When employing a beta diversity metric, use the weighted
#'        version. Default: \code{TRUE}.
#'        
#' @param rank   What rank of taxa to display, for example \code{"Phylum"} or 
#'        \code{"Genus"}. Run \code{taxa_ranks()} to see all options for a 
#'        given BIOM object. The default, \code{NULL}, selects the lowest
#'        level.
#'        
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{5}.
#'
#' @param p.top   Only display taxa with the most significant differences in
#'        abundance. If \code{p.top} is >= 1, then the \code{p.top} most 
#'        significant taxa are displayed. If \code{p.top} is less than one, all
#'        taxa with an adjusted p-value <= \code{p.top} are displayed.
#'        Recommended to be used in combination with the \code{taxa} parameter
#'        to set a lower bound on the mean abundance of considered taxa.
#'        Default: \code{Inf}.
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of p-values.
#'        Run \code{p.adjust.methods} for a list of available options.
#'        Default: \code{"fdr"}.
#'
#' @param perms   Number of random permutations to use for estimating statistical
#'        significance. Default: \code{1000}.
#'        
#' @param ...   Parameters for underlying functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, \code{dot.size = 2} or \code{d.size = 2} ensures only the 
#'        dotplot layer has its size set to \code{2}.
#'        
#' @return A \code{ggplot2} plot. The computed data points and statistics will 
#'         be attached as \code{attr(p, 'data')} and \code{attr(p, 'stats')}, 
#'         respectively.
#' 
#' 
#' @export
#' @seealso \code{\link{stats_table}}
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50) 
#'     bdiv_biplot(biom, color.by="Body Site", rank="Genus")
#'     
bdiv_biplot <- function (
  biom, metric = "Bray-Curtis", ord = "PCoA", layers = NULL, 
  color.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
  weighted = TRUE, rank = NULL, taxa = 5, p.top = Inf, p.adj = "fdr", 
  perms = 1000, ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("bdiv_biplot", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  history <- attr(biom, 'history')
  history %<>% c(sprintf("bdiv_biplot(%s)", as.args(params, fun = bdiv_biplot)))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    metric %<>% validate_arg(biom, n = 1, 'metric', 'bdiv', tree = tree)
    ord    %<>% validate_arg(biom, n = 1, 'ord')
    rank   %<>% validate_arg(biom, n = 1, 'rank', default = tail(c('OTU', taxa_ranks(biom)), 1))
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
  if (nsamples(biom) < 4)
    stop("At least four samples are needed for an ordination.")
  
  
  #________________________________________________________
  # Convert user's `layers` spec to layer names.
  #________________________________________________________
  layer_names <- with(params, {
    
    layerlist <- c( 'p' = "point", 'n' = "name",    's' = "spider", 
                    'n' = "name",  'd' = "density", 't' = "taxon",
                    'a' = "arrow", 'e' = "ellipse", 'm' = "mean" )
    
    default <- if (is_null(color.by)) "p" else c("p", "c", "e")
    if (!is_null(rank)) default %<>% c("t", "m")
    
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
    biom     = biom, 
    dist     = params[['metric']],
    ord      = params[['ord']],
    weighted = params[['weighted']],
    md       = TRUE,
    split.by = params[['facet.by']],
    stat.by  = names(params[['color.by']]),
    perms    = params[['perms']],
    rank     = params[['rank']],
    taxa     = params[['taxa']],
    p.adj    = params[['p.adj']],
    p.top    = params[['p.top']] )
  
  ggdata %<>% rename_cols(
    ".axis.1" = ".x", 
    ".axis.2" = ".y",
    ".sample" = ".label" )
  
  
  
  
  #________________________________________________________
  # Initialize the `layers` object
  #________________________________________________________
  layers <- list()
  
  attr(layers, 'biom')     <- biom
  attr(layers, 'data')     <- ggdata
  attr(layers, 'params')   <- params
  attr(layers, 'function') <- bdiv_biplot
  attr(layers, 'xcol')     <- ".x"
  attr(layers, 'ycol')     <- ".y"
  attr(layers, 'xmode')    <- "numeric"
  attr(layers, 'ymode')    <- "numeric"
  
  if (!is_null(params[['color.by']]))   layer_names %<>% c("color")
  if (!is_null(params[['shape.by']]))   layer_names %<>% c("shape")
  if (!is_null(params[['pattern.by']])) layer_names %<>% c("pattern")
  
  initLayer(layer_names)
  # layers <- metadata_layers(layers)
  
  remove("layer_names")
  
  
  
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
  # Statistics table to show the user
  #________________________________________________________
  attr(p, 'stats') <- list('groupwise' = attr(ggdata, "stats_tbl", exact = TRUE))
  
  
  set_cache_value(cache_file, p)
  return(p)
}




#________________________________________________________
# Edit ggdata and biplot for ggplot compatibility
#________________________________________________________
ordination_biplot <- function (layers) {
  
  if (!any(c("mean", "arrow", "taxon") %in% names(layers)))
    return (layers)
  
  ggdata <- attr(layers, 'data',   exact = TRUE)
  biplot <- attr(ggdata, 'biplot', exact = TRUE)
  params <- attr(layers, 'params', exact = TRUE)
  
  if (is_null(biplot))           return (layers)
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
  facet_cols <- c(params[['facet.by']], '.weight', '.dist', '.ord', '.rank')
  
  biplot <- data.frame(
    check.names = FALSE,
    biplot[,facet_cols,drop=FALSE],
    '.label' = biplot[['.taxa']],
    '.x'     = biplot[['.axis.1']],
    '.y'     = biplot[['.axis.2']],
    '.xend'  = biplot[['.ori.1']],
    '.yend'  = biplot[['.ori.2']],
    '.p.val' = biplot[['.p.val']],
    '.adj.p' = biplot[['.adj.p']],
    '.size'  = biplot[['.abundance']] )
    
  
  attr(ggdata, 'biplot') <- biplot
  attr(layers, 'data')   <- ggdata
  
  return (layers)
}




#________________________________________________________
# Facetting logic.
#________________________________________________________
ordination_facets <- function (layers) {
  
  ggdata <- attr(layers, 'data',      exact = TRUE)
  params <- attr(layers, 'params',    exact = TRUE)
  stats  <- attr(ggdata, 'stats_tbl', exact = TRUE)
  biplot <- attr(ggdata, 'biplot',    exact = TRUE)
  ranks  <- params[['rank']]
  
  
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
  for (i in c('.weight', '.dist', '.ord')) {
    df[[i]] %<>% as.character()
    if (length(unique(df[[i]])) > 1)
      var_metrics %<>% c(i)
  }
  
  if (length(var_metrics) == 0) {
    
    # Plot Title = "Weighted UniFrac PCoA (Phylum)"
    metric <- paste(collapse = " ", df[1, c('.weight', '.dist', '.ord')])
    if (length(ranks) == 1) metric <- sprintf("%s (%s)", metric, ranks)
    setLayer("labs", title = metric)
    
    
  } else if (length(var_metrics) == 1) {
    
    # Plot Title = "Weighted UniFrac (Phylum)"
    metric <- paste(collapse = " ", df[1, setdiff(c('.weight', '.dist', '.ord'), var_metrics)])
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
      yes  = sprintf("%s<br>%s %s %s", df[['.facet']], df[['.weight']], df[['.dist']], df[['.ord']]), 
      no   = sprintf("**%s %s %s**",                   df[['.weight']], df[['.dist']], df[['.ord']]) )
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
  if (!is_null(stats)) {
    
    # stats are agnostic of .ord and .rank
    df <- plyr::join(df, stats, by = c('.weight', '.dist', params[['facet.by']]))
    
    stats_text <- sprintf(
      fmt = "*p* = %s; *R<sup>2</sup>* = %s; *F* = %s",
      format(df[['.p.val']],  digits=3), 
      format(df[['.r.sqr']],  digits=3), 
      format(df[['.f.stat']], digits=3) )
    
    if (length(unique(stats_text)) == 1) {
      setLayer("labs",  subtitle      = stats_text[[1]])
      setLayer("theme", plot.subtitle = element_markdown())
      
    } else {
      df[['.facet']] <- bool_switch(
        test = hasName(df, ".facet"),
        yes  = sprintf("%s<br>%s", df[['.facet']], stats_text), 
        no   = stats_text )
    }
    
    setLayer("labs", caption = local({
      sprintf("Statistics computed with Adonis (%i permutations).", params[['perms']]) }))
    setLayer("theme", plot.caption = element_text(face = "italic"))
  }
  
  
  #________________________________________________________
  # Copy .facet column from df to biplot
  #________________________________________________________
  if (!is_null(biplot) && hasName(df, ".facet"))
    attr(df, 'biplot') <- local({
      join_cols <- c('.weight', '.dist', '.ord', '.rank', params[['facet.by']])
      join_cols <- intersect(join_cols, intersect(colnames(df), colnames(biplot)))
      
      biplot <- plyr::join(
        x     = biplot, 
        y     = df[,c(join_cols, '.facet'),drop=FALSE], 
        by    = join_cols, 
        match = "first" )
      
      for (i in join_cols)
        biplot[[i]] %<>% factor(levels = levels(ggdata[[i]]))
      biplot[['.facet']] %<>% factor(levels = unique(df[['.facet']]))
      
      return (biplot)
    })
  
  
  #________________________________________________________
  # Convert back to factors with original levels ordering.
  #________________________________________________________
  for (i in c('.weight', '.dist', '.ord', '.rank', params[['facet.by']]))
    if (hasName(df, i)) df[[i]] %<>% factor(levels = levels(ggdata[[i]]))
  
  
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


