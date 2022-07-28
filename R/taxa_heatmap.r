#' Display taxa abundances as a heatmap.
#' 
#' @name taxa_heatmap
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
#' @param taxa   Which taxa to give separate rows. An integer value will show
#'        the top n most abundant taxa. A value 0 <= n < 1 will show all taxa 
#'        with that mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{6}.
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
#'     taxa_heatmap(biom, rank="Phylum", color.by="Body Site")
#'     
taxa_heatmap <- function (
    biom, rank = NULL, taxa = 6,
    color.by = NULL, order.by = NULL, 
    colors = NULL, orders = NULL, ...) {
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  rank  <- bool_switch(
    test = is.null(rank), 
    yes  = tail(c('OTU', taxa_ranks(biom)), 1), 
    no   = as.vector(validate_metrics(biom, rank, "rank", multi=TRUE)) )
  
  if (!is.null(order.by))
    order.by <- as.vector(validate_metrics(biom, order.by, "meta"))
  
  
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
      do.call(taxa_heatmap, args)
    })
    
    p <- patchwork::wrap_plots(plots, ncol = 1)
    
    history <- sprintf("taxa_heatmap(%s)", as.args(params, fun = taxa_heatmap))
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
  # Subset before calculating relative abundances
  #________________________________________________________
  if (!is.null(params[['order.by']]))
    metadata(biom) <- subset_by_params(
      df     = metadata(biom),
      params = params[c('order.by', 'orders')] )
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (nsamples(biom) < 1)
    stop("At least one sample is needed for a taxa heatmap.")
  
  
  
  #________________________________________________________
  # Matrix with samples in columns and taxa in rows.
  #________________________________________________________
  mtx <- t(taxa_rollup(
    biom = biom, 
    rank = params[['rank']], 
    taxa = params[['taxa']] ))
  
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
  
  if (is.null(params[['colors']])) {
    args[['colors']] <- "khroma::bilbao"
    names(args[['colors']]) <- sprintf(".%s Abundance", params[['rank']])
  } else {
    # Inject some defaults without interfering with user's specs.
    args[['colors']] <- structure(
      .Data      = params[['colors']],
      grid_name  = sprintf("%s Abundance", params[['rank']]),
      grid_color = "khroma::bilbao" )
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
  history <- sprintf("taxa_heatmap(%s)", as.args(params, fun = taxa_heatmap))
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return(p)
  
}

