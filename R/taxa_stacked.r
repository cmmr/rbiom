#' Display taxa abundances as a stacked bar graph.
#' 
#' @inherit documentation_default
#' @inherit documentation_plot_return return
#' 
#' @family taxa_abundance
#' @family visualization
#' 
#' @param colors,patterns   A character vector of colors or patterns to use in
#'        the graph. A named character vector can be used to map taxon names to 
#'        specific colors or patterns. Set to `TRUE` to auto-select colors
#'        or patterns, or to `FALSE` to disable per-taxa colors or patterns. 
#'        Default: \code{colors=TRUE, patterns=FALSE}.
#' 
#' @param label.by,order.by   What metadata column to use for labeling and/or
#'        sorting the samples across the x-axis. Set \code{label.by='.sample'} 
#'        to display sample names. When \code{order.by=NULL}, samples are 
#'        arranged based on \code{dist} and \code{clust}, below.
#'        Default: \code{label.by=NULL, order.by=NULL}.
#' 
#' @param dist,clust   Distance ([stats::dist()]) and clustering 
#'        ([stats::hclust()]) methods to use for automatically arranging 
#'        samples along the x-axis to put samples with similar composition near 
#'        one another. Default: \code{dist="euclidean", clust="complete"}.
#' 
#' @param ...   Parameters for underlying functions. Prefixing parameter names 
#'        with a layer name ensures that a particular parameter is passed to, 
#'        and only to, that layer.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     taxa_stacked(biom, rank="Phylum")
#'     
#'     taxa_stacked(biom, rank = "genus", facet.by = "body site")
#'     
taxa_stacked <- function (
    biom, rank = -1, taxa = 6, colors = TRUE, patterns = FALSE,
    label.by = NULL, order.by = NULL, facet.by = NULL, 
    dist = "euclidean", clust = "complete", other = TRUE, 
    unc = "singly", lineage = FALSE, xlab.angle = 90, ...) {
  
  biom <- as_rbiom(biom)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env(..., .dots = TRUE)
  cache_file <- get_cache_file('taxa_stacked', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  params <- list2env(params)
  
  
  #________________________________________________________
  # Handle multiple ranks with recursive subcalls.
  #________________________________________________________
  validate_rank(max = Inf, env = params)
  
  if (length(params$rank) > 1) {
    
    ranks <- params$rank
    
    plots <- sapply(ranks, simplify = FALSE, function (rank) {
      params[['rank']]       <- rank
      params[['labs.title']] <- rank
      do.call(taxa_stacked, fun_params(taxa_stacked, params))
    })
    
    cmd <- current_cmd('taxa_stacked')
    p   <- plot_wrap(plots, cmd, ncol = 1)
    
    set_cache_value(cache_file, p)
    return (p)
  }
  
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  with(params, {
    
    validate_dist()
    validate_clust()
    
    validate_biom_field('label.by', null_ok = TRUE)
    validate_biom_field('order.by', null_ok = TRUE, max = Inf)
    validate_biom_field('facet.by', null_ok = TRUE, max = Inf, col_type = "cat")
    
    na.omit(biom, clone = FALSE, fields = c(label.by, order.by, facet.by))
    
    if (biom$n_samples < 1)
      stop("At least one sample is needed for a stacked bar plot.")
  })
  
  
  
  
  #________________________________________________________
  # Compute taxa abundance values.
  #________________________________________________________
  with(params, {
    
    .ggdata <- local({
      
      
      # The matrix will be needed for similarity clustering.
      #________________________________________________________
      mtx <- taxa_matrix(
        biom    = biom, 
        rank    = rank, 
        taxa    = taxa, 
        lineage = lineage,
        sparse  = FALSE, 
        unc     = unc, 
        other   = other )
      
      
      
      # Pivot to long-form and add metadata.
      #________________________________________________________
      tbl <- tibble(
          '.rank'      = rank,
          '.sample'    = colnames(mtx)[col(mtx)],
          '.taxa'      = rownames(mtx)[row(mtx)] %>% factor(rownames(mtx)),
          '.abundance' = as.numeric(mtx) ) %>%
        left_join(biom$metadata, by = '.sample') %>%
        add_class('rbiom_tbl')
      
      
      
      # Factor levels on .sample control x-axis ordering.
      #________________________________________________________
      tbl[['.sample']] %<>% factor(levels = local({
        
        dm   <- dist(x = t(mtx), method = dist)
        hc   <- hclust(d = dm, method = clust)
        lvls <- hc$labels[hc$order]
        
        if (!is.null(order.by)) {
          
          sns <- biom$samples
          md  <- biom$metadata[,order.by]
          md[['.sort']] <- match(sns, lvls)
          
          lvls <- sns[base::rank(do.call(order, unname(as.list(md))))]
        }
        
        return (lvls)
        
      }))
      
      
      
      return (tbl)
    })
    
    
    .xcol  <- '.sample'
    .ycol  <- '.abundance'
    .xmode <- 'factor'
    
    stat.by    <- ".taxa"
    color.by   <- ".taxa"
    pattern.by <- ".taxa"
  })
  
  
  
  #________________________________________________________
  # Provenance-tracked ggplot2 functions.
  #________________________________________________________
  .element_blank <- P('ggplot2::element_blank')
  .element_text  <- P('ggplot2::element_text')
  .unit          <- P('grid::unit')
  .label_number  <- P('scales::label_number')
  .cut_si        <- P('scales::cut_si')
  
  
  
  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  init_layers(params, do_init = "stack")
  
  
  
  #________________________________________________________
  # Display a label on the x-axis other than sample name.
  #________________________________________________________
  if (is.null(params$label.by)) {
    
    set_layer(params, 'xaxis', labels = NULL, breaks = NULL)
    
  } else if (!eq(params$label.by, ".sample")) {
    
    # Convert `.sample` to "<int>-<label>" format for easy parsing.
    with(params, {
      levels(.ggdata[['.sample']]) <- local({
        lvls <- levels(.ggdata[['.sample']])
        labs <- pull(biom, label.by)
        paste(sep="-", seq_along(lvls), labs[lvls])
      })
    })
    
    set_layer(params, 'xaxis', labels = ~ sub("^\\d+\\-", "", .))
  }
  
  
  #________________________________________________________
  # Caption about labeling and ordering
  #________________________________________________________
  set_layer(
    params = params, 
    layer  = 'labs', 
    'x'       = "Sample",
    'caption' = with(params, glue(
    
      if ((is_null(label.by) || eq(label.by, '.sample')) && is_null(order.by)) {
        "Arranged by {dist} distance and '{clust}' clustering."
        
      } else if ((is_null(label.by) || eq(label.by, '.sample'))) {
        "Ordered by {order.by}."
        
      } else if (is_null(order.by)) {
        "Labeled by {label.by} and arranged by\n{dist} distance and\n'{clust}' clustering."
        
      } else {
        ifelse(
          test = eq(label.by, order.by), 
          yes  = "Labeled and ordered by {label.by}.", 
          no   = "Labeled by {label.by} and ordered by {order.by}." )
      }
    )))
  set_layer(params, 'theme', plot.caption = .element_text(size = 9, face = "italic"))
  
  
  #________________________________________________________
  # If rarefied, display relative abundance in %
  #________________________________________________________
  depth <- unique(sample_sums(params$biom))
  if (length(depth) == 1) {
    
    params$.ggdata[[params$.ycol]] <- params$.ggdata[[params$.ycol]] / depth
    
    set_layer(
      params = params, 
      layer  = 'labs', 
      'y' = "Relative Abundance" )
    
    set_layer(
      params = params, 
      layer  = 'yaxis',
      'expand' = c(0,0),
      'limits' = c(0,1),
      'labels' = scales::percent )
    
    set_layer(
      params = params, 
      layer  = 'theme', 
      'panel.grid'  = .element_blank(),
      'plot.margin' = .unit(c(1,1,1,1), "lines") )
    
    
  } else {
    
    set_layer(
      params = params, 
      layer  = 'labs',  
      'y' = "Raw Abundance" )
    
    set_layer(
      params = params, 
      layer  = 'yaxis', 
      'expand' = c(0, 0, 0.02, 0),
      'breaks' = with(params, base::pretty(.ggdata[[.ycol]])),
      'labels' = .label_number(scale_cut = .cut_si("")) )
    
    set_layer(
      params = params, 
      layer  = 'theme',
      'panel.grid.major.x' = .element_blank(),
      'panel.grid.minor.y' = .element_blank() )
  }
  
  
  
  #________________________________________________________
  # aes() parameters
  #________________________________________________________
  set_layer(
    params = params, 
    layer  = 'stack', 
    'mapping|x' = params$.xcol, 
    'mapping|y' = params$.ycol )
  
  
  
  
  #________________________________________________________
  # Non-aes parameters
  #________________________________________________________
  legend_title <- ifelse(eq(params$labs.title, params$rank), "", params$rank)
  if (!(is.null(params$colors) || isFALSE(params$colors))) {
    set_layer(params, 'labs', color = legend_title)
    set_layer(params, 'labs', fill  = legend_title)
    
  } else if (!(is.null(params$patterns) || isFALSE(params$patterns))) {
    set_layer(params, 'labs', fill  = legend_title)
  }
  remove("legend_title")
  
  if (isTRUE(params$xlab.angle == 30)) {
    set_layer(params, 'theme', 'axis.text.x' = .element_text(angle=-30, vjust=1, hjust=0))
  } else {
    set_layer(params, 'theme', 'axis.text.x' = .element_text(angle=-90, vjust=0.3, hjust=0))
  }
    
  
  
  
  #________________________________________________________
  # Create the plot and add each layer with its arguments.
  # Also attaches a human-readable version of the plot command.
  #________________________________________________________
  p <- params %>%
    plot_facets() %>%
    plot_build()
  
  
  attr(p, 'cmd') <- current_cmd('taxa_stacked')
  set_cache_value(cache_file, p)
  
  return (p)
}


