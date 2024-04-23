#' Display beta diversities in an all vs all grid.
#' 
#' @inherit documentation_heatmap
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family visualization
#' 
#' @param grid   Color palette name, or a list with entries for \code{label}, 
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or 
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: `"devon"`
#'                 
#' @param color.by   Add annotation tracks for these metadata column(s). 
#'        See "Annotation Tracks" section below for details.
#'        Default: `NULL`
#'                 
#' @param order.by   Which metadata column(s) to use for ordering the samples 
#'        across the x and y axes. Overrides any \code{clust} argument.
#'        See "Ordering and Limiting" section below for details.
#'        Default: `NULL`
#'                 
#' @param limit.by   Metadata definition(s) to use for sample subsetting prior
#'        to calculations. 
#'        See "Ordering and Limiting" section below for details.
#'        Default: `NULL`
#'        
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#'        For example, \code{labs.subtitle = "Plot subtitle"}.
#'         
#'         
#' @section Annotation Tracks:
#' 
#' Metadata can be displayed as colored tracks above the heatmap. Common use 
#' cases are provided below, with more thorough documentation available at 
#' https://cmmr.github.io/rbiom .
#' 
#' \preformatted{## Categorical ----------------------------
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
#' \preformatted{## Categorical ----------------------------
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
#' \preformatted{## Categorical ----------------------------
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
#'     library(rbiom)
#'     
#'     # Keep and rarefy the 10 most deeply sequenced samples.
#'     hmp10 <- rarefy(hmp50, n = 10)
#'     
#'     bdiv_heatmap(hmp10, color.by=c("Body Site", "Age"))
#'     
#'     bdiv_heatmap(hmp10, bdiv="uni", weighted=c(T,F), color.by="sex")
#'     
bdiv_heatmap <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, grid = "devon",
    color.by = NULL, order.by = NULL, limit.by = NULL, 
    label = TRUE, label_size = NULL, rescale = "none", clust = "complete", 
    trees = TRUE, asp = 1, tree_height = 10, track_height = 10, 
    legend = "right", title = TRUE, xlab.angle = "auto", ...) {
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("bdiv_heatmap(%s)", as.args(params, 2, bdiv_heatmap))
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_heatmap', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
   with(params, {
    
    validate_bdiv(max = Inf)
    
    validate_bool("trees")
    validate_bool("weighted", max = Inf)
    
    validate_meta_aes('color.by', null_ok = TRUE, max = Inf, aes = "color")
    validate_meta_aes('order.by', null_ok = TRUE, max = Inf)
    validate_meta_aes('limit.by', null_ok = TRUE, max = Inf)
    
    if (!is_list(grid)) grid <- list(label = "Dissimilarity", colors = grid)
    if (length(order.by) > 0) clust <- FALSE
  })
  
  
  #________________________________________________________
  # Handle multiple metrics/weightings with recursion.
  #________________________________________________________
  if (length(params$bdiv) > 1 || length(params$weighted) > 1) {
    
    plots <- list()
    
    for (d in params$bdiv)
      for (w in params$weighted)
        plots[[length(plots) + 1]] <- local({
          
          args <- fun_params(bdiv_heatmap, params)
          
          args[['bdiv']]     <- d
          args[['weighted']] <- w
          
          if (isTRUE(args[['title']]))
            args[['title']] <- paste(ifelse(w, "Weighted", "Unweighted"), d)
          
          do.call(bdiv_heatmap, args)
        })
    
    p <- patchwork::wrap_plots(plots)
    
    attr(p, 'data') <- lapply(plots, attr, which = 'data', exact = TRUE)
    
    attr(p, 'code') <- paste(collapse = "\n\n", local({
      
      cmds <- sapply(seq_along(plots), function (i) {
        sub(
          x           = attr(plots[[i]], 'code', exact = TRUE), 
          pattern     = "ggplot(data)", 
          replacement = sprintf("p%i <- ggplot(data[[%i]])", i, i),
          fixed       = TRUE )
      })
      c(cmds, sprintf("patchwork::wrap_plots(%s)", paste0(collapse = ", ", "p", seq_along(plots))))
      
    })) %>% add_class('rbiom_code')
    
    return (p)
  }
  
  
  #________________________________________________________
  # Subset biom by requested metadata and aes.
  #________________________________________________________
  sync_metadata(params)
  biom <- params$biom
  
  
  #________________________________________________________
  # Sanity Check
  #________________________________________________________
  if (biom$n_samples < 1)
    stop("At least one sample is needed for a bdiv heatmap.")
  
  
  #________________________________________________________
  # Replace title=TRUE with a default title string.
  #________________________________________________________
  with(params, {
    
    if (isTRUE(title))
      title <- ifelse(
        test = isTRUE(weighted), 
        yes  = paste("Weighted",   bdiv, "Distance"), 
        no   = paste("Unweighted", bdiv, "Distance") )
    
    if (!(is.null(title) || (is_scalar_character(title) && !is.na(title))))
      cli_abort("title must be TRUE, NULL, or a character string, not {.type {title}}: {title}.")
  })
    
  
  
  
  #________________________________________________________
  # Matrix of samples x samples.
  #________________________________________________________
  dm  <- bdiv_distmat(
    biom     = biom, 
    bdiv     = params$bdiv, 
    weighted = params$weighted, 
    tree     = params$tree )
  mtx <- as.matrix(dm)
  
  
  
  #________________________________________________________
  # Arguments to pass on to plot_heatmap
  #________________________________________________________
  args <- fun_params(plot_heatmap, params)
  args[['mtx']]  <- mtx
  args[['dist']] <- dm
  
  
  for (md_col in names(params$color.by))
    params$color.by[[md_col]] %<>% within({
      if (exists("values", inherits = FALSE))
        colors <- values
      values <- pull(biom, md_col)
    })
  args[['tracks']] <- params$color.by
  
  
  
  fig <- do.call(plot_heatmap, args)
  
  attr(fig, 'cmd') <- cmd
  set_cache_value(cache_file, fig)
  
  return (fig)
}

