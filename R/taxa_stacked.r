#' Display taxa abundances as a stacked bar graph.
#' 
#' @name taxa_stacked
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#'        
#' @param rank   What rank(s) of taxa to display, for example \code{"Phylum"} 
#'        or \code{"Genus"}. Run \code{taxa_ranks()} to see all options for a 
#'        given BIOM object. The default, \code{NULL}, selects the lowest
#'        level.
#'        
#' @param taxa   Which taxa to give separate colors. An integer value will show
#'        the top n most abundant taxa. A value 0 <= n < 1 will show all taxa 
#'        with that mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{6}.
#'                 
#' @param colors,patterns   A character vector of colors or patterns to use in
#'        the graph. A named character vector can be used to map taxon names to 
#'        specific colors or patterns. Set to \code{TRUE} to auto-select colors
#'        or patterns, or to \code{FALSE} to disable per-taxa colors
#'        or patterns. Default: \code{colors=TRUE, patterns=FALSE}.
#'                 
#' @param facet.by   Split into multiple plots based on this metadata 
#'        column(s). Default: \code{NULL}.
#'                 
#' @param label.by,order.by   What metadata column to use for labeling and/or
#'        sorting the samples across the x-axis. When \code{order.by=NULL},
#'        samples are arranged based on \code{dist} and \code{clust}, below.
#'        Default: \code{label.by=NULL, order.by=NULL}.
#'                 
#' @param facets,labels,orders   Subset to just these values in the 
#'        \code{facet.by}, \code{label.by}, and/or \code{order.by} metadata
#'        columns. Default: \code{facets=NULL, labels=NULL, orders=NULL}.
#'                 
#' @param dist,clust   Distance (\link{dist}) and clustering (\link{hclust})
#'        methods to use for automatically arranging samples along the x-axis 
#'        to group samples with similar composition near one another.
#'        Default: \code{dist="euclidean", clust="complete"}.
#'                 
#' @param other   Add an 'Other' taxa that ensures all bars sum to 100%. 
#'        Default: \code{TRUE}.
#'        
#' @param xlab.angle   How to rotate the tick labels on the x-axis.
#'        Options are \code{90} (the default) or \code{30}.
#'        
#' @param ...   Parameters for underlying functions. Prefixing parameter names 
#'        with a layer name ensures that a particular parameter is passed to, 
#'        and only to, that layer.
#'        
#' @return A \code{ggplot2} plot. The computed data points will be attached as 
#'         \code{attr(, 'data')}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     taxa_stacked(biom, rank="Phylum")
#'     
taxa_stacked <- function (
    biom, rank = NULL, taxa = 6, 
    colors = TRUE, patterns = FALSE,
    facet.by = NULL, label.by = NULL, order.by = NULL, 
    facets = NULL, labels = NULL, orders = NULL, 
    dist = "euclidean", clust = "complete",
    other = TRUE, xlab.angle = 90, ...) {
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  dist  <- as.vector(validate_metrics(NULL, dist,  "dist"))
  clust <- as.vector(validate_metrics(NULL, clust, "clust"))
  rank  <- bool_switch(
    test = is.null(rank), 
    yes  = tail(c('OTU', taxa_ranks(biom)), 1), 
    no   = as.vector(validate_metrics(biom, rank, "rank", multi=TRUE)) )
  
  if (!is.null(facet.by)) facet.by <- as.vector(validate_metrics(biom, facet.by, "meta", multi=TRUE))
  if (!is.null(label.by)) label.by <- as.vector(validate_metrics(biom, label.by, "meta"))
  if (!is.null(order.by)) order.by <- as.vector(validate_metrics(biom, order.by, "meta"))
  
  
  #________________________________________________________
  # Collect all parameters into a list
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  remove(list = setdiff(ls(), c("params", "biom")))
  
  
  #________________________________________________________
  # Handle multiple ranks with recursive subcalls.
  #________________________________________________________
  if (length(params[['rank']]) > 1) {
    
    ranks <- params[['rank']]
    
    plots <- sapply(ranks, simplify = FALSE, function (rank) {
      args                 <- params
      args[['rank']]       <- rank
      args[['labs.title']] <- rank
      do.call(taxa_stacked, args)
    })
    
    p <- patchwork::wrap_plots(plots, ncol = 1)
     
    history <- sprintf("taxa_stacked(%s)", as.args(params, fun = taxa_stacked))
    attr(p, 'history') <- c(attr(biom, 'history'), history)
    attr(p, 'data')    <- lapply(plots, attr, which = 'data', exact = TRUE)
    attr(p, 'cmd')     <- paste(collapse = "\n\n", local({
      cmds <- sapply(seq_along(ranks), function (i) {
        sub(
          x           = attr(plots[[i]], 'cmd', exact = TRUE), 
          pattern     = "ggplot(data)", 
          replacement = sprintf("p%i <- ggplot(data[[%s]])", i, glue::single_quote(ranks[[i]])),
          fixed       = TRUE )
      })
      c(cmds, sprintf("patchwork::wrap_plots(%s, ncol = 1)", paste0(collapse = ", ", "p", seq_along(ranks))))
    }))
    
    return (p)
  }
  
  
  #________________________________________________________
  # Store the cmd history before script-modding the params.
  #________________________________________________________
  history <- sprintf("taxa_stacked(%s)", as.args(params, fun = taxa_stacked))
  
  
  #________________________________________________________
  # Inject color.by and pattern.by arguments as needed.
  #________________________________________________________
  taxa <- structure(".taxa", display=FALSE)
  if (!isFALSE(params[['colors']]))   params[['color.by']]   <- taxa
  if (isTRUE(params[['colors']]))     params[['colors']]     <- NULL
  if (!isFALSE(params[['patterns']])) params[['pattern.by']] <- taxa
  if (isTRUE(params[['patterns']]))   params[['patterns']]   <- NULL
  remove("taxa")
  
  
  #________________________________________________________
  # Subset before calculating relative abundances
  #________________________________________________________
  metadata(biom) %<>% subset_by_params(params)
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (nsamples(biom) < 1)
    stop("At least one sample is needed for a stacked bar plot.")
  
  
  
  #________________________________________________________
  # Limit to top n taxa, min abundance, or custom list
  #________________________________________________________
  ggdata <- local({
    
    taxa <- params[['taxa']]
    mat  <- t(taxa_matrix(biom, rank = params[['rank']]))
  
    if (is.numeric(taxa)) {
      rel <- sort(rowMeans(t(t(mat) / colSums(mat))), decreasing = TRUE)
      if (taxa >= 1) { taxa <- head(names(rel), taxa) 
      } else         { taxa <- names(rel)[rel >= taxa] }
      
    } else {
      taxa <- intersect(as.character(taxa), rownames(mat))
    }
    
    if (length(taxa) == 0)
      stop("No taxa match the criteria: ", capture.output(str(taxa)))
    
    #________________________________________________________
    # Sum the remaining taxa into an "Other" category
    #________________________________________________________
    if (length(taxa) == nrow(mat) || isFALSE(params[['other']])) {
      mat <- mat[taxa,,drop=FALSE]
      
    } else {
      other <- matrix(
        data = colSums(mat[!rownames(mat) %in% taxa,,drop=FALSE]), 
        nrow = 1, dimnames = list("Other", colnames(mat)))
      mat <- rbind(mat[taxa,,drop=FALSE], other)
    }
    
      
    #________________________________________________________
    # Cluster samples by composition similarity or order.by
    #________________________________________________________
    sorted <- local({
      if (is.null(params[['order.by']])) {
        dm <- dist(x = t(mat), method = params[['dist']])
        hc <- hclust(d = dm, method = params[['clust']])
        setNames(seq_along(hc$labels), hc$labels[hc$order])
      } else {
        rank(metadata(biom, params[['order.by']]))
      }
    })
    
    
    #________________________________________________________
    # Convert counts matrix to a data.frame
    #________________________________________________________
    data.frame(
      check.names = FALSE,
      stringsAsFactors = FALSE,
      '.sample' = colnames(mat)[col(mat)],
      '.rank'   = params[['rank']],
      '.taxa'   = rownames(mat)[row(mat)],
      '.value'  = as.numeric(mat),
      '.sort'   = sorted[colnames(mat)[col(mat)]] )
  })
  
  
  #________________________________________________________
  # Convert taxa names to a factor
  #________________________________________________________
  vals <- sort(unique(ggdata[['.taxa']]))
  taxa <- params[['taxa']]
  taxa <- if (is.character(taxa)) intersect(taxa, vals) else vals
  if (isTRUE(params[['other']]) && "Other" %in% vals)
    taxa <- c(setdiff(taxa, "Other"), "Other")
  ggdata[['.taxa']] %<>% factor(levels = taxa)
  
  remove("vals", "taxa")
  
  
  
  #________________________________________________________
  # Order samples on x-axis by similarity or order.by
  #________________________________________________________
  lvls <- ggdata[!duplicated(ggdata[['.sample']]), , drop=FALSE]
  lvls <- lvls[order(lvls[['.sort']]), '.sample', drop=TRUE]
  ggdata[['.sample']] %<>% factor(levels = lvls)
  remove("lvls")
  
  
  
  #________________________________________________________
  # Initialize the `layers` object
  #________________________________________________________
  layers <- list()
  
  attr(layers, 'biom')     <- biom
  attr(layers, 'data')     <- ggdata
  attr(layers, 'params')   <- params
  attr(layers, 'function') <- taxa_stacked
  attr(layers, 'xcol')     <- ".sample"
  attr(layers, 'ycol')     <- ".value"
  attr(layers, 'xmode')    <- "factor"
  attr(layers, 'ymode')    <- "numeric"
  
  initLayer(c('ggplot', 'stack', 'labs', 'theme', 'theme_bw'))
  layers <- metadata_layers(layers)
  layers[['color']] <- NULL
  remove("ggdata")
  
  
  #________________________________________________________
  # Add in facet.by columns
  #________________________________________________________
  if (!is.null(params[['facet.by']])) {
    ggdata <- attr(layers, 'data', exact = TRUE)
    
    ids <- as.character(ggdata[['.sample']])
    for (i in params[['facet.by']])
      ggdata[[i]] <- metadata(biom, i)[ids]
    
    attr(layers, 'data') <- ggdata
    remove("ggdata", "ids", "i")
    
    setLayer(layer = "facet", scales = "free_x")
  }
  
  
  
  #________________________________________________________
  # Display a label on the x-axis other than sample name.
  #________________________________________________________
  if (!is.null(params[['label.by']])) {
    
    # Convert `.sample` to "<int>-<label>" format for easy parsing.
    ggdata <- attr(layers, 'data', exact = TRUE)
    levels(ggdata[['.sample']]) <- local({
      levels <- levels(ggdata[['.sample']])
      labels <- metadata(biom, params[['label.by']])
      paste(sep="-", seq_along(levels), labels[levels])
    })
    attr(layers, 'data') <- ggdata
    remove("ggdata")
    
    setLayer(layer = "xaxis", 'labels' = ~ sub("^\\d+\\-", "", .))
  }
  
  
  #________________________________________________________
  # Caption about labeling and ordering
  #________________________________________________________
  setLayer(
    layer = "labs", 
    'x'       = "Sample",
    'caption' = with(params, glue::glue(
    
      if (is.null(label.by) && is.null(order.by)) {
        "Samples are ordered according to '{clust}' clustering based on {dist} distance."
        
      } else if (is.null(label.by)) {
        "Samples are ordered by {order.by})*"
        
      } else if (is.null(order.by)) {
        "Samples are labeled by {label.by} and ordered according\nto '{clust}' clustering based on {dist} distance."
        
      } else {
        ifelse(
          test = identical(label.by, order.by), 
          yes  = "Samples are labeled and ordered by {label.by}.", 
          no   = "Samples are labeled by {label.by} and ordered by {order.by}." )
      }
    )))
  setLayer(layer = "theme", 'plot.caption' = element_text(face = "italic", color = "gray"))
  
  
  #________________________________________________________
  # Control the values displayed on the y axis
  #________________________________________________________
  # browser()
  if (is_rarefied(biom)) {
    ggdata <- attr(layers, 'data', exact = TRUE)
    ggdata[['.value']] <- ggdata[['.value']] / attr(biom, 'rarefaction', exact = TRUE)
    
    setLayer(
      layer = "yaxis",
      'expand' = c(0,0),
      'limits' = c(0,1),
      'labels' = scales::percent )
    
    setLayer(
      layer = "theme", 
      'panel.grid'  = element_blank(),
      'plot.margin' = as.cmd(unit(c(1,1,1,1), "lines")) )
    
    setLayer(
      layer = "labs", 
      'y' = "Relative Abundance" )
    
    attr(layers, 'data') <- ggdata
    remove("ggdata")
    
  } else {
    
    ggdata <- attr(layers, 'data', exact = TRUE)
    breaks <- base::pretty(ggdata[['.value']])
    labels <- siunit(breaks)
    
    setLayer(
      layer = "labs",  
      'y' = "Raw Abundance" )
    setLayer(
      layer = "yaxis", 
      'expand' = c(0, 0, 0.02, 0),
      'breaks' = breaks,
      'labels' = labels )
    setLayer(
      layer = "theme",
      'panel.grid.major.x' = element_blank(),
      'panel.grid.minor.y' = element_blank() )
    
    remove("ggdata", "breaks", "labels")
  }
  
  
  
  #________________________________________________________
  # aes() parameters
  #________________________________________________________
  setLayer(layer = "stack", 'mapping|x' = ".sample", 'mapping|y' = ".value")
  
  
  #________________________________________________________
  # Pattern-specific aes / non-aes mappings
  #________________________________________________________
  layer <- "stack"
  if (hasName(layers, "pattern")) {
    setLayer(
      'mapping|pattern_type' = ".taxa", 
      'pattern'              = "magick", 
      'pattern_res'          = 120,
      'fill'                 = "white",
      'color'                = "black" )
    
    if (isFALSE(params[['colors']])) { setLayer('pattern_fill'         = "black")
    } else                           { setLayer('mapping|pattern_fill' = ".taxa") }
    
  } else {
    setLayer('mapping|fill' = ".taxa")
  }
  remove("layer")
  
  
  
  
  #________________________________________________________
  # Non-aes parameters
  #________________________________________________________
  legend_title <- ifelse(identical(params[['labs.title']], params[['rank']]), "", params[['rank']])
  if (hasName(layers, "pattern")) {
                                  setLayer(layer = "labs", pattern_type  = legend_title)
    if (hasName(layers, "fill"))  setLayer(layer = "labs", pattern_fill  = legend_title)
    if (hasName(layers, "color")) setLayer(layer = "labs", pattern_color = legend_title)
    
  } else {
    if (hasName(layers, "fill"))  setLayer(layer = "labs", fill  = legend_title)
    if (hasName(layers, "color")) setLayer(layer = "labs", color = legend_title)
  }
  remove("legend_title")
  
  if (isTRUE(params[['xlab.angle']] == 30)) {
    setLayer(layer = "theme", 'axis.text.x' = element_text(angle=-30, vjust=1, hjust=0))
  } else {
    setLayer(layer = "theme", 'axis.text.x' = element_text(angle=-90, vjust=0.3, hjust=0))
  }
    
  
  
  
  #________________________________________________________
  # Create the plot and add each layer with its arguments.
  # Also attaches a human-readable version of the plot command.
  #________________________________________________________
  p <- layers %>%
    plot_facets() %>%
    plot_build()
  
  
  #________________________________________________________
  # Attach history of biom modifications and this call
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return(p)
  
}


