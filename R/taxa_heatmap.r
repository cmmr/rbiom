#' Display taxa abundances as a heatmap.
#' 
#' @inherit documentation_heatmap
#' @inherit documentation_default
#' @inherit bdiv_heatmap sections
#' 
#' @family taxa_abundance
#' @family visualization
#' 
#' @param grid   Color palette name, or a list as expected [plot_heatmap()].
#'        Default: `"bilbao"`
#' 
#' @param tracks   A character vector of metadata fields to display as tracks 
#'        at the top of the plot. Or, a list as expected by the `tracks` 
#'        argument of [plot_heatmap()]. Default: `NULL`
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Keep and rarefy the 10 most deeply sequenced samples.
#'     hmp10 <- rarefy(hmp50, n = 10)
#'     
#'     taxa_heatmap(hmp10, rank = "Phylum", tracks = "Body Site")
#'     
#'     taxa_heatmap(hmp10, rank = "Genus", tracks = c("sex", "bo"))
#'     
#'     taxa_heatmap(hmp10, rank = "Phylum", tracks = list(
#'       'Sex'       = list(colors = c(m = "#0000FF", f = "violetred")), 
#'       'Body Site' = list(colors = "muted", label = "Source") ))
#'     
taxa_heatmap <- function (
    biom, rank = -1, taxa = 6, tracks = NULL, grid = "bilbao",
    other = FALSE, unc = "singly", lineage = FALSE, 
    label = TRUE, label_size = NULL, rescale = "none", trees = TRUE,
    clust = "complete", dist = "euclidean", 
    asp = 1, tree_height = 10, track_height = 10, 
    legend = "right", title = TRUE, xlab.angle = "auto", ...) {
  
  biom <- as_rbiom(biom)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env(...)
  cache_file <- get_cache_file('taxa_heatmap', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  params <- list2env(params)
  
  
  #________________________________________________________
  # Handle multiple ranks with recursive sub-calls.
  #________________________________________________________
  validate_rank(max = Inf, env = params)
  
  if (length(params$rank) > 1) {
    
    ranks <- params$rank
    
    plots <- sapply(ranks, simplify = FALSE, function (rank) {
      
      args <- fun_params(taxa_heatmap, params)
      args[['rank']] <- rank
      if (is.null(params$title))
        args[['title']] <- TRUE
      
      do.call(taxa_heatmap, args)
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
    
    attr(p, 'cmd') <- current_cmd('bdiv_heatmap')
    
    set_cache_value(cache_file, p)
    return (p)
  }
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  with(params, {
    
    validate_taxa()
    # validate_biom_field('order.by', null_ok = TRUE, max = Inf)
    validate_var_choices('rescale', choices = c("none", "cols", "rows"))
    
    if (isFALSE(title)) title <- NULL
    if (isTRUE(title))  title <- paste(rank, "Abundance")
    validate_var_length('title', title, max = 1, null_ok = TRUE)
    
    if (biom$n_samples < 1)
      cli_abort("At least one sample is needed for a heatmap.")
    
  })
  
  
  
  #________________________________________________________
  # Default values
  #________________________________________________________
  with(params, {
    
    if (!is_list(grid)) grid <- list(colors = grid)
    if (!is_scalar_character(grid$label)) {
      grid$label <- if (rescale != "none") { "Rescaled Abundance"
      } else if (is.null(biom$depth))      { "Unrarefied Reads"
      } else if (biom$depth == 1)          { "Relative Abundance"
      } else                               { "Rarefied Reads" }
    }
    
    if (length(clust)    < 1) clust <- "complete"
    # if (length(order.by) > 0) clust <- c(clust[[1]], NA)
  })
  
  
  
  #________________________________________________________
  # Convert `tracks` into a named list of lists.
  #________________________________________________________
  biom_tracks(params)
  
  
  
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
    
    remove("biom", "rank", "taxa", "lineage", "unc", "other")
  })
  
  
  
  #________________________________________________________
  # Actual plotting is handled by plot_heatmap()
  #________________________________________________________
  p <- do.call(plot_heatmap, as.list(params))
  
  attr(p, 'cmd') <- current_cmd('taxa_heatmap')
  set_cache_value(cache_file, p)
  
  return (p)
}

