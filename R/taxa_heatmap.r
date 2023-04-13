#' Display taxa abundances as a heatmap.
#' 
#' @name taxa_heatmap
#' 
#' @inherit plot_heatmap params return
#' @inherit bdiv_heatmap params sections
#' 
#' @family plotting
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
#' @param grid   Color palette name, or a list with entries for \code{label}, 
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or 
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: \code{list(label = "{rank} Abundance", colors = "bilbao")}.
#'        
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#'        For example, \code{labs.title = "Plot Title"}.
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
    grid = list(label = "{rank} Abundance", colors = "bilbao"),
    color.by = NULL, order.by = NULL, limit.by = NULL, 
    label = TRUE, label_size = NULL, rescale = "none", trees = TRUE,
    clust = "complete", dist = "euclidean", 
    tree_height = NULL, track_height = NULL, ratio=1, 
    legend = "right", xlab.angle = "auto", ...) {
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- attr(biom, 'history')
  history %<>% c(sprintf("taxa_heatmap(%s)", as.args(params, fun = taxa_heatmap)))
  remove(list = setdiff(ls(), c("params", "history")))
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    stopifnot(is_scalar_atomic(taxa))
    rank %<>% validate_arg(biom, 'rank', n = c(1,Inf), default = tail(c('OTU', taxa_ranks(biom)), 1))
    
    if (!is_list(grid)) grid <- list(label = "{rank} Abundance", colors = grid)
    if (length(clust)    < 1) clust <- "complete"
    if (length(order.by) > 0) clust <- c(clust[[1]], NA)
  })
  
  
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
    
    attr(p, 'history') <- history
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
    
    remove("ranks", "plots")
    
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
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (nsamples(biom) < 1)
    stop("At least one sample is needed for a taxa heatmap.")
  
  
  
  #________________________________________________________
  # Matrix with samples in columns and taxa in rows.
  #________________________________________________________
  mtx <- t(taxa_matrix(
    biom = biom, 
    rank = params[['rank']], 
    taxa = params[['taxa']] ))
  
  
  
  #________________________________________________________
  # Arguments to pass on to plot_heatmap
  #________________________________________________________
  excl <- setdiff(formalArgs(taxa_heatmap), formalArgs(plot_heatmap))
  excl <- c(excl, names(params)[startsWith(names(params), '.')])
  args <- params[setdiff(names(params), excl)]
  
  args[['mtx']] <- mtx
  args[['grid']][['label']] %<>% sub("{rank}", params[['rank']], ., fixed = TRUE)
  
  for (md_col in names(params[['color.by']]))
    params[['color.by']][[md_col]] %<>% within({
      colors <- values
      values <- metadata(biom, md_col)
    })
  args[['tracks']] <- params[['color.by']]
  
  
  p <- do.call(plot_heatmap, args)
  
  attr(p, 'history') <- history
  
  
  return(p)
}

