#' Display beta diversities in an all vs all grid.
#' 
#' @name bdiv_heatmap
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
#' @param color.by   Add annotation tracks for these metadata column(s). 
#'        Default: \code{NULL}.
#'                 
#' @param order.by   Which metadata columns to use for ordering the samples 
#'        across the x-axis. The default setting, \code{order.by=NULL}, uses 
#'        sample names as labels and arranges samples by similarity. See 
#'        \link{plot_heatmap} for details.
#'                 
#' @param colors   A named \code{list()} with names matching \code{color.by}, 
#'        each specifying a \bold{palette} as described in \link{plot_heatmap}.
#'        For example: \code{colors = list(Sex = c(Male="blue", Female="red"))}.
#'        Default: \code{NULL}.
#'                 
#' @param orders   Subset to just these values in the \code{order.by} metadata 
#'        column. Default: \code{NULL}.
#'        
#' @param xlab.angle   How to rotate the tick labels on the x-axis. 
#'        \bold{'auto'} (the default), automatically selects a rotation value. 
#'        \bold{0}, \bold{30}, and \bold{90} sets the angle to horizontal, 
#'        angled, and vertical, respectively.
#'        
#' @param tree   A phylogenetic tree for use in calculating UniFrac distance.
#'        The default, \code{NULL}, will use the BIOM object's tree.
#'                 
#' @param ...   Additional arguments to pass on to \link{plot_heatmap}.
#'        
#' @return A \code{ggplot2} plot. The computed data points will be attached as 
#'         \code{attr(, 'data')}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50 %>% rarefy() %>% select(1:10)
#'     bdiv_heatmap(biom, color.by="Body Site")
#'     
bdiv_heatmap <- function (
    biom, metric = "Bray-Curtis", weighted = TRUE,
    color.by = NULL, order.by = NULL, 
    colors = NULL, orders = NULL,
    xlab.angle = "auto", tree = NULL, ...) {
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  
  if (!is.null(order.by))
    order.by <- as.vector(validate_metrics(biom, order.by, "meta", multi=TRUE))
  
  if (!is.null(color.by))
    color.by <- as.vector(validate_metrics(biom, color.by, "meta", multi=TRUE))
  
  
  #________________________________________________________
  # Collect all parameters into a list
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  remove(list = setdiff(ls(), c("params", "biom")))
  
  
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
      c(cmds, sprintf("patchwork::wrap_plots(%s)", paste0(collapse = ", ", "p", seq_along(ranks))))
    }))
    
    return (p)
  }
  
  
  #________________________________________________________
  # Subset
  #________________________________________________________
  if (!is.null(params[['order.by']]) || !is.null(params[['color.by']]))
    metadata(biom) <- subset_by_params(
      df     = metadata(biom),
      params = params[c('order.by', 'orders', 'color.by', 'colors')] )
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (nsamples(biom) < 1)
    stop("At least one sample is needed for a bdiv heatmap.")
  
  
  
  #________________________________________________________
  # Matrix of samples x samples.
  #________________________________________________________
  mtx <- as.matrix(bdiv_dist(
    biom     = biom, 
    method   = params[['metric']], 
    weighted = params[['weighted']],
    tree     = params[['tree']] ))
  
  # supports multi-column ordering
  if (!is.null(params[['order.by']]))
    mtx <- mtx[,do.call(order, lapply(
      X   = params[['order.by']], 
      FUN = function (i) {
        as.vector(metadata(biom, i)[colnames(mtx)]) }))]
  
  
  
  #________________________________________________________
  # Arguments to pass on to plot_heatmap
  #________________________________________________________
  args <- list(...)
  
  args[['mtx']] <- mtx
  
  within(args, {
    labs.title %||=% ifelse(
      test = isTRUE(params[['weighted']]), 
      yes  = paste0("Weighted\n",   params[['metric']], "\nDistance"), 
      no   = paste0("Unweighted\n", params[['metric']], "\nDistance") )})
  
  if (is.null(params[['colors']])) {
    args[['colors']] <- c(.Distance = "-khroma::bilbao")
  } else {
    # Inject some defaults without interfering with user's specs.
    args[['colors']] <- structure(
      .Data      = params[['colors']],
      grid_name  = "Distance",
      grid_color = "khroma::bilbao",
      grid_rev   = TRUE )
  }
  
  if (length(params[['order.by']]) > 0 && is.null(args[['clust']]))
    args[['clust']] <- c("complete", NA)
  
  if (!is.null(params[['color.by']])) {
    df <- metadata(biom)
    df <- df[colnames(mtx), params[['color.by']], drop=FALSE]
    args[['top_tracks']] <- df
    
    remove("df")
  }
  
  p <- do.call(plot_heatmap, args)
  
  
  #________________________________________________________
  # Attach history of biom modifications and this call
  #________________________________________________________
  history <- sprintf("bdiv_heatmap(%s)", as.args(params, fun = bdiv_heatmap))
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return(p)
  
}

