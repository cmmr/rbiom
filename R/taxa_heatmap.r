#' Display taxa abundances as a heatmap.
#' 
#' @inherit documentation_heatmap
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' @family visualization
#' 
#' @param grid   Color palette name, or a list with entries for \code{label}, 
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or 
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: \code{list(label = "{rank} Abundance", colors = "bilbao")}.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Keep and rarefy the 10 most deeply sequenced samples.
#'     hmp10 <- rarefy(hmp50, n = 10)
#'     
#'     taxa_heatmap(hmp10, rank="Phylum", color.by="Body Site")
#'     
taxa_heatmap <- function (
    biom, rank = -1, taxa = 6,
    grid = list(label = "{rank} Abundance", colors = "bilbao"),
    color.by = NULL, order.by = NULL, limit.by = NULL, 
    other = FALSE, unc = "singly", lineage = FALSE, 
    label = TRUE, label_size = NULL, rescale = "none", trees = TRUE,
    clust = "complete", dist = "euclidean", 
    tree_height = NULL, track_height = NULL, ratio=1, 
    legend = "right", xlab.angle = "auto", ...) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment(), ...)
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
      do.call(taxa_heatmap, fun_params(taxa_heatmap, params))
    })
    
    p <- patchwork::wrap_plots(plots, ncol = 1) %>% 
      add_class('rbiom_plot')
    
    attr(p, 'data') <- lapply(plots, attr, which = 'data', exact = TRUE)
    
    attr(p, 'code') <- paste(collapse = "\n\n", local({
      cmds <- sapply(seq_along(ranks), function (i) {
        sub(
          x           = plots[[i]][['ggcmd']], 
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
  # Default values
  #________________________________________________________
  with(params, {
    
    if (!is_list(grid))                   grid       <- list(colors = grid)
    if (!is_scalar_character(grid$label)) grid$label <- "{rank} Abundance"
    grid$label %<>% sub("{rank}", rank, ., fixed = TRUE)
    
    if (length(clust)    < 1) clust <- "complete"
    if (length(order.by) > 0) clust <- c(clust[[1]], NA)
  })
  
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    
    validate_taxa()
    
    validate_meta_aes('color.by', null_ok = TRUE, max = Inf)
    validate_meta_aes('order.by', null_ok = TRUE, max = Inf)
    validate_meta_aes('limit.by', null_ok = TRUE, max = Inf)
    
    sync_metadata()
    
    if (biom$n_samples < 1)
      stop("At least one sample is needed for a heatmap.")
  })
  
  
  
  #________________________________________________________
  # Matrix with samples in columns and taxa in rows.
  #________________________________________________________
  with(params, {
    
    mtx <- taxa_matrix(
      biom    = biom, 
      rank    = rank, 
      taxa    = taxa,
      lineage = lineage, 
      sparse  = FALSE, 
      unc     = unc, 
      other   = other )
    
  })
  
  
  
  #________________________________________________________
  # Convert aes `color.by` spec to `tracks` format.
  #________________________________________________________
  for (md_col in names(params$color.by))
    params$color.by[[md_col]] %<>% within({
      if (exists('values', inherits = FALSE))
        colors <- values
      values <- pull(params$biom, md_col)
    })
  params[['tracks']] <- params$color.by
  
  
  
  #________________________________________________________
  # Arguments to pass on to plot_heatmap
  #________________________________________________________
  args <- intersect(formalArgs(plot_heatmap), env_names(params))
  args <- c(mget(args, envir = params), params$.dots)
  p    <- do.call(plot_heatmap, args)
  
  
  set_cache_value(cache_file, p)
  return (p)
}

