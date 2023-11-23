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
#'        specific colors or patterns. Set to \code{TRUE} to auto-select colors
#'        or patterns, or to \code{FALSE} to disable per-taxa colors
#'        or patterns. Default: \code{colors=TRUE, patterns=FALSE}.
#'        
#' @param facet.by,limit.by   Metadata columns to 
#'        use for data partitioning. Default: \code{NULL}
#'                 
#' @param label.by,order.by   What metadata column to use for labeling and/or
#'        sorting the samples across the x-axis. When \code{order.by=NULL},
#'        samples are arranged based on \code{dist} and \code{clust}, below.
#'        Default: \code{NULL}.
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
#'     biom <- sample_rarefy(hmp50)
#'     taxa_barplot(biom, rank="Phylum")
#'     
taxa_barplot <- function (
    biom, rank = -1, taxa = 6, colors = TRUE, patterns = FALSE,
    label.by = NULL, order.by = NULL, facet.by = NULL, limit.by = NULL, 
    dist = "euclidean", clust = "complete", other = TRUE, 
    unc = "singly", lineage = FALSE, xlab.angle = 90, ...) {
  
  validate_biom(clone = FALSE)
  
  params  <- eval_envir(environment(), ...)
  history <- append_history('fig ', params)
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Handle multiple ranks with recursive subcalls.
  #________________________________________________________
  validate_rank(max = Inf, env = params)
  
  if (length(params$rank) > 1) {
    
    ranks <- params$rank
    
    plots <- sapply(ranks, simplify = FALSE, function (rank) {
      params[['rank']]       <- rank
      params[['labs.title']] <- rank
      do.call(taxa_barplot, fun_params(taxa_barplot, params))
    })
    
    p <- patchwork::wrap_plots(plots, ncol = 1) %>% 
      add_class('rbiom_plot')
    
    attr(p, 'history') <- history
    attr(p, 'data')    <- lapply(plots, attr, which = 'data', exact = TRUE)
    
    attr(p, 'code') <- paste(collapse = "\n\n", local({
      cmds <- sapply(seq_along(ranks), function (i) {
        sub(
          x           = attr(plots[[i]], 'code', exact = TRUE), 
          pattern     = "ggplot(data)", 
          replacement = sprintf("p%i <- ggplot(data[[%s]])", i, single_quote(ranks[[i]])),
          fixed       = TRUE )
      })
      c(cmds, sprintf("patchwork::wrap_plots(%s, ncol = 1)", paste0(collapse = ", ", "p", seq_along(ranks))))
      
    })) %>% add_class('rbiom_code')
    
    set_cache_value(cache_file, p)
    return (p)
  }
  
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  with(params, {
    
    validate_dist()
    validate_clust()
    
    validate_meta_aes('label.by', null_ok = TRUE)
    validate_meta_aes('order.by', null_ok = TRUE, max = Inf)
    validate_meta_aes('facet.by', null_ok = TRUE, max = Inf, col_type = "cat")
    validate_meta_aes('limit.by', null_ok = TRUE, max = Inf)
    
    sync_metadata()
    
    if (n_samples(biom) < 1)
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
        left_join(sample_metadata(biom), by = '.sample') %>%
        add_class('rbiom_tbl')
      
      
      
      # Factor levels on .sample control x-axis ordering.
      #________________________________________________________
      tbl[['.sample']] %<>% factor(levels = local({
        
        dm   <- dist(x = t(mtx), method = dist)
        hc   <- hclust(d = dm, method = clust)
        lvls <- hc$labels[hc$order]
        
        if (!is.null(order.by)) {
          
          sns <- sample_names(biom)
          md  <- sample_metadata(biom)[,order.by]
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
    
  })
  
  
  
  #________________________________________________________
  # Inject color.by and pattern.by arguments as needed.
  #________________________________________________________
  
  taxa <- levels(params$.ggdata[['.taxa']])
  n    <- length(taxa)
  
  if (!isFALSE(values <- params$colors)) {
    if (!is_null(names(values))) { values <- params$colors[taxa]
    } else                       { values <- get_n_colors(n, values) }
    params$color.by <- list('.taxa' = list('values' = values))
    attr(params$color.by, 'display') <- FALSE
  }
  
  if (!isFALSE(values <- params$patterns)) {
    if (!is_null(names(values))) { values <- params$patterns[taxa]
    } else                       { values <- get_n_patterns(n, values) }
    params$pattern.by <- list('.taxa' = list('values' = values))
    attr(params$pattern.by, 'display') <- FALSE
  }
  
  remove("taxa", "n")
  
  
  
  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  do_init <- "stack"
  if (!is_null(params$color.by))   do_init %<>% c("fill")
  if (!is_null(params$pattern.by)) do_init %<>% c("pattern")
  if (!is_null(params$facet.by))   do_init %<>% c("facet")
  
  init_layers(params, do_init = do_init)
  
  
  
  #________________________________________________________
  # Display a label on the x-axis other than sample name.
  #________________________________________________________
  if (!is_null(params$label.by)) {
    
    # Convert `.sample` to "<int>-<label>" format for easy parsing.
    with(params, {
      levels(.ggdata[['.sample']]) <- local({
        lvls <- levels(.ggdata[['.sample']])
        labs <- sample_metadata(biom, label.by)
        paste(sep="-", seq_along(lvls), labs[lvls])
      })
    })
    
    set_layer(params, 'xaxis', 'labels' = ~ sub("^\\d+\\-", "", .))
  }
  
  
  #________________________________________________________
  # Caption about labeling and ordering
  #________________________________________________________
  set_layer(
    params = params, 
    layer  = 'labs', 
    'x'        = "Sample",
    'subtitle' = with(params, glue(
    
      if (is_null(label.by) && is_null(order.by)) {
        "Samples are ordered according to '{clust}' clustering based on {dist} distance."
        
      } else if (is_null(label.by)) {
        "Samples are ordered by {order.by}."
        
      } else if (is_null(order.by)) {
        "Samples are labeled by {label.by} and ordered according\nto '{clust}' clustering based on {dist} distance."
        
      } else {
        ifelse(
          test = eq(label.by, order.by), 
          yes  = "Samples are labeled and ordered by {label.by}.", 
          no   = "Samples are labeled by {label.by} and ordered by {order.by}." )
      }
    )))
  set_layer(params, 'theme', 'plot.subtitle' = element_text(size = 9, lineheight = 1.2))
  
  
  #________________________________________________________
  # Control the values displayed on the y axis
  #________________________________________________________
  
  if (is_rarefied(params$biom)) {
    
    with(params, {
      .ggdata[[.ycol]] <- .ggdata[[.ycol]] / rare_depth(biom)
    })
    
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
      'panel.grid'  = element_blank(),
      'plot.margin' = as.cmd(unit(c(1,1,1,1), "lines")) )
    
    
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
      'labels' = si_units )
    
    set_layer(
      params = params, 
      layer  = 'theme',
      'panel.grid.major.x' = element_blank(),
      'panel.grid.minor.y' = element_blank() )
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
  # Pattern-specific aes / non-aes mappings
  #________________________________________________________
  if (has_layer(params, 'pattern')) {
    
    set_layer(
      params = params, 
      layer  = 'stack',
      'mapping|pattern_type' = ".taxa", 
      'pattern'              = "magick", 
      'pattern_res'          = 120,
      'fill'                 = "white",
      'color'                = "black" )
    
    if (isFALSE(params$colors)) { set_layer(params, 'stack', 'pattern_fill'         = "black")
    } else                      { set_layer(params, 'stack', 'mapping|pattern_fill' = ".taxa") }
    
  } else {
    set_layer(params, 'stack', 'mapping|fill' = ".taxa")
  }
  
  
  
  
  #________________________________________________________
  # Non-aes parameters
  #________________________________________________________
  legend_title <- ifelse(eq(params$labs.title, params$rank), "", params$rank)
  if (has_layer(params, 'pattern')) {
                                    set_layer(params, 'labs', pattern_type  = legend_title)
    if (has_layer(params, 'fill'))  set_layer(params, 'labs', pattern_fill  = legend_title)
    if (has_layer(params, 'color')) set_layer(params, 'labs', pattern_color = legend_title)
    
  } else {
    if (has_layer(params, 'fill'))  set_layer(params, 'labs', fill  = legend_title)
    if (has_layer(params, 'color')) set_layer(params, 'labs', color = legend_title)
  }
  remove("legend_title")
  
  if (isTRUE(params$xlab.angle == 30)) {
    set_layer(params, 'theme', 'axis.text.x' = element_text(angle=-30, vjust=1, hjust=0))
  } else {
    set_layer(params, 'theme', 'axis.text.x' = element_text(angle=-90, vjust=0.3, hjust=0))
  }
    
  
  
  
  #________________________________________________________
  # Create the plot and add each layer with its arguments.
  # Also attaches a human-readable version of the plot command.
  #________________________________________________________
  p <- params %>%
    plot_facets() %>%
    plot_build()
  
  attr(p, 'history') <- history
  
  
  set_cache_value(cache_file, p)
  return (p)
}

