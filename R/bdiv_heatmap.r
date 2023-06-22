#' Display beta diversities in an all vs all grid.
#' 
#' @name bdiv_heatmap
#' 
#' @inherit plot_heatmap params return
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param metric   Beta diversity metric(s) to use. Options are 
#'        \code{"Manhattan"}, \code{"Euclidean"}, \code{"Bray-Curtis"}, 
#'        \code{"Jaccard"}, and \code{"UniFrac"}. UniFrac requires a 
#'        phylogenetic tree. Default: \code{"Bray-Curtis"}.
#'        
#' @param weighted   Should the beta diversity metric be weighted by 
#'        abundances? Default: \code{TRUE}.
#' 
#' @param grid   Color palette name, or a list with entries for \code{label}, 
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or 
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: \code{list(label = "Distance", colors = "-bilbao")}.
#'                 
#' @param color.by   Add annotation tracks for these metadata column(s). 
#'        See "Annotation Tracks" section below for details.
#'        Default: \code{NULL}
#'                 
#' @param order.by   Which metadata column(s) to use for ordering the samples 
#'        across the x and y axes. Overrides any \code{clust} argument.
#'        See "Ordering and Limiting" section below for details.
#'        Default: \code{NULL}
#'                 
#' @param limit.by   Metadata definition(s) to use for sample subsetting prior
#'        to calculations. 
#'        See "Ordering and Limiting" section below for details.
#'        Default: \code{NULL}
#'        
#' @param tree   A phylogenetic tree for use in calculating UniFrac distance.
#'        The default, \code{NULL}, will use the BIOM object's tree.
#'        
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#'        For example, \code{labs.title = "Plot Title"}.
#'        
#' @return A \code{ggplot2} plot. The computed data points will be attached as 
#'         \code{attr(, 'data')}.
#'         
#'         
#' @section Annotation Tracks:
#' 
#' Metadata can be displayed as colored tracks above the heatmap. Common use 
#' cases are provided below, with more thorough documentation available at 
#' https://cmmr.github.io/rbiom .
#' 
#' \preformatted{  ## Categorical ----------------------------
#' color.by = "Body Site"
#' color.by = list('Body Site' = "bright")
#' color.by = list('Body Site' = c("Stool", "Saliva"), 'colors' = "bright")
#' color.by = list('Body Site' = c('Stool' = "blue", 'Saliva' = "green"))
#' 
#' ## Numeric --------------------------------
#' color.by = "Age"
#' color.by = list('Age' = "reds")
#' color.by = list('Age' = c(20,NA), 'colors' = "reds") # at least 20 years old
#' color.by = list('Age' = c(20,40)) # between 20 and 40 years old (inclusive)
#' 
#' ## Multiple Tracks ------------------------
#' color.by = c("Body Site", "Age")
#' color.by = list('Body Site' = "bright", 'Age' = "reds")
#' color.by = list(
#'   'Body Site' = c('Stool' = "blue", 'Saliva' = "green"),
#'   'Age'       = list(range = c(20,40), 'colors' = "reds") )
#' }
#' 
#' The following entries in the track definitions are understood: 
#' 
#' \itemize{
#'   \item{\code{colors} - }{ A pre-defined palette name or custom set of colors to map to. }
#'   \item{\code{range} - }{ The c(min,max) to use for scale values. }
#'   \item{\code{label} - }{ Label for this track. Defaults to the name of this list element. }
#'   \item{\code{side} - }{ Options are \code{"top"} (default) or \code{"left"}. }
#'   \item{\code{na.color} - }{ The color to use for \code{NA} values. }
#'   \item{\code{bins} - }{ Bin a gradient into this many bins/steps. }
#'   \item{\code{guide} - }{ A list of arguments for guide_colorbar() or guide_legend(). }
#' }
#' 
#' All built-in color palettes are colorblind-friendly.
#' 
#' Categorical palette names: \code{"okabe"}, \code{"carto"}, \code{"r4"}, 
#' \code{"polychrome"}, \code{"tol"}, \code{"bright"}, \code{"light"}, 
#' \code{"muted"}, \code{"vibrant"}, \code{"tableau"}, \code{"classic"}, 
#' \code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, and \code{"fishy"}.
#' 
#' Numeric palette names: \code{"reds"}, \code{"oranges"}, \code{"greens"}, 
#' \code{"purples"}, \code{"grays"}, \code{"acton"}, \code{"bamako"}, 
#' \code{"batlow"}, \code{"bilbao"}, \code{"buda"}, \code{"davos"}, 
#' \code{"devon"}, \code{"grayC"}, \code{"hawaii"}, \code{"imola"}, 
#' \code{"lajolla"}, \code{"lapaz"}, \code{"nuuk"}, \code{"oslo"}, 
#' \code{"tokyo"}, \code{"turku"}, \code{"bam"}, \code{"berlin"}, 
#' \code{"broc"}, \code{"cork"}, \code{"lisbon"}, \code{"roma"}, 
#' \code{"tofino"}, \code{"vanimo"}, and \code{"vik"}.
#'         
#'         
#' @section Ordering and Limiting:
#' 
#' \bold{\code{order.by}} controls which metadata column(s) are used to arrange
#' samples on the plot. It also enables subsetting to a particular set or 
#' range of values. Prefix a column name with \code{-} to arrange values in 
#' descending order rather than ascending.
#' 
#' \preformatted{  ## Categorical ----------------------------
#' order.by = "Body Site"
#' order.by = list('Body Site' = c("Stool", "Saliva"))
#' 
#' ## Numeric --------------------------------
#' order.by = "-Age"
#' order.by = list('Age'  = c(20,NA)) # at least 20 years old
#' order.by = list('-Age' = c(20,40)) # between 20 and 40 years old (inclusive)
#' 
#' ## Multiple / Mixed -----------------------
#' order.by = c("-Body Site", "Age")
#' order.by = list("Body Site", '-Age' = c(20,40))
#' }
#' 
#' 
#' \bold{\code{limit.by}} is used to specify a subset of samples without any
#' side-effects on aesthetics. It is especially useful for limiting the data
#' to a single categorical metadata value. Unlike the other *.by parameters,
#' \code{limit.by} must always be a named \code{list()}.
#' 
#' \preformatted{  ## Categorical ----------------------------
#' limit.by = list('Sex' = "Male")
#' 
#' ## Numeric --------------------------------
#' limit.by = list('Age' = c(20,NA)) # at least 20 years old
#' limit.by = list('Age' = c(20,40)) # between 20 and 40 years old (inclusive)
#' 
#' ## Multiple / Mixed -----------------------
#' limit.by = list(
#'   'Sex'       = "Male", 
#'   'Body Site' = c("Stool", "Saliva")
#'   'Age'       = c(20,40) )
#' }
#' 
#' 
#' 
#' 
#' @export
#' @examples
#'   library(rbiom) 
#'   
#'   biom <- hmp50 %>% rarefy() %>% select(1:10)
#'   bdiv_heatmap(biom, color.by="Body Site")
#'     
bdiv_heatmap <- function (
    biom, metric = "Bray-Curtis", weighted = TRUE,
    grid = list(label = "Distance", colors = "-bilbao"),
    color.by = NULL, order.by = NULL, limit.by = NULL, 
    label = TRUE, label_size = NULL, rescale = "none", trees = TRUE,
    clust = "complete", dist = "euclidean", 
    tree_height  = NULL, track_height = NULL, ratio=1, 
    legend = "right", xlab.angle = "auto", ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("bdiv_heatmap", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  history <- attr(biom, 'history')
  history %<>% c(sprintf("bdiv_heatmap(%s)", as.args(params, fun = bdiv_heatmap)))
  remove(list = setdiff(ls(), c("params", "history", "cache_file")))
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    stopifnot(is_scalar_logical(weighted))
    metric %<>% validate_arg(biom, 'metric', 'bdiv', n = c(1,Inf))
    
    if (!is_list(grid)) grid <- list(label = "Distance", colors = grid)
    if (length(order.by) > 0) clust <- FALSE
  })
  
  
  #________________________________________________________
  # Handle multiple ranks with recursive subcalls.
  #________________________________________________________
  if (length(params[['metric']]) > 1 || length(params[['weighted']]) > 1) {
    
    plots <- list()
    
    for (m in params[['metric']])
      for (w in params[['weighted']])
        plots[[length(plots) + 1]] <- local({
          args                 <- params
          args[['metric']]     <- m
          args[['weighted']]   <- w
          args[['labs.title']] <- paste(ifelse(w, "Weighted", "Unweighted"), m)
          do.call(bdiv_heatmap, args)
        })
    
    p <- patchwork::wrap_plots(plots)
    
    history <- sprintf("bdiv_heatmap(%s)", as.args(params, fun = bdiv_heatmap))
    attr(p, 'history') <- history
    attr(p, 'data')    <- lapply(plots, attr, which = 'data', exact = TRUE)
    attr(p, 'cmd')     <- paste(collapse = "\n\n", local({
      cmds <- sapply(seq_along(ranks), function (i) {
        sub(
          x           = attr(plots[[i]], 'cmd', exact = TRUE), 
          pattern     = "ggplot(data)", 
          replacement = sprintf("p%i <- ggplot(data[[%s]])", i, single_quote(ranks[[i]])),
          fixed       = TRUE )
      })
      c(cmds, sprintf("patchwork::wrap_plots(%s)", paste0(collapse = ", ", "p", seq_along(ranks))))
    }))
    
    return (p)
  }
  
  
  #________________________________________________________
  # Subset biom by requested metadata and aes.
  #________________________________________________________
  params %<>% metadata_params(contraints = list(
    color.by = list(),
    order.by = list(),
    limit.by = list() ))
  
  biom <- params[['biom']]
  
  
  # #________________________________________________________
  # # Subset
  # #________________________________________________________
  # if (!is_null(params[['order.by']]) || !is_null(params[['color.by']]))
  #   metadata(biom) <- subset_by_params(
  #     df     = metadata(biom),
  #     params = params[c('order.by', 'orders', 'color.by', 'colors')] )
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (nsamples(biom) < 1)
    stop("At least one sample is needed for a bdiv heatmap.")
  
  
  
  #________________________________________________________
  # Matrix of samples x samples.
  #________________________________________________________
  mtx <- as.matrix(bdiv_distmat(
    biom     = biom, 
    method   = params[['metric']], 
    weighted = params[['weighted']] ))
  
  
  
  #________________________________________________________
  # Arguments to pass on to plot_heatmap
  #________________________________________________________
  excl <- setdiff(formalArgs(bdiv_heatmap), formalArgs(plot_heatmap))
  excl <- c(excl, names(params)[startsWith(names(params), '.')])
  args <- params[setdiff(names(params), excl)]
  args[['mtx']] <- mtx
  
  within(args, {
    labs.title %||=% ifelse(
      test = isTRUE(params[['weighted']]), 
      yes  = paste0("Weighted\n",   params[['metric']], "\nDistance"), 
      no   = paste0("Unweighted\n", params[['metric']], "\nDistance") )})
  
  for (md_col in names(params[['color.by']]))
    params[['color.by']][[md_col]] %<>% within({
      colors <- values
      values <- metadata(biom, md_col)
    })
  args[['tracks']] <- params[['color.by']]
  
  
  p <- do.call(plot_heatmap, args)
  
  attr(p, 'history') <- history
  
  
  set_cache_value(cache_file, p)
  return (p)
}

