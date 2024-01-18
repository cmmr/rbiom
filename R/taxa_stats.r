

#' Test taxa abundances for significant differences.
#' 
#' @inherit documentation_test.ifelse
#' @inherit documentation_model.lm
#' @inherit documentation_default
#' @inherit documentation_stats_return return
#' 
#' @family taxa_abundance
#' @family stats_tables
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     taxa_stats(biom, stat.by = "Body Site", rank = "Family")
#'     
#'     taxa_stats(biom, stat.by = "Body Site", regr = "Age", rank = "Family")

taxa_stats <- function (
    biom, stat.by = NULL, regr = NULL, rank = -1, taxa = 6, 
    test = ifelse(is.null(regr), "means", "trends"), 
    model = "lm", level = 0.95, 
    trans = ifelse(is.null(regr), "none", "rank"), 
    lineage = FALSE, unc = "singly", other = FALSE,
    split.by = NULL, p.adj = "fdr" ) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment())
  remove(list = intersect(env_names(params), ls()))
  
  with(params, {
    md       <- setdiff(unique(c(stat.by, regr, split.by)), ".taxa")
    split.by <- setdiff(c(split.by, ".taxa"), c(regr, stat.by))
  })
  
  params$df <- do.call(taxa_table,  fun_params(taxa_table,  params))
  stats     <- do.call(stats_table, fun_params(stats_table, params))
  
  
  return (stats)
}




#' Apply `p.top` constraint to data about to be plotted.
#' 
#' Used by boxplot_stats() and corrplot_stats().
#' 
#' @noRd
#' @keywords internal
#' 

apply_p.top <- function (params) {
  
  p.top  <- params$p.top
  ggdata <- params$.ggdata
  stats  <- params$.plot_attrs[['stats']]
  vline  <- attr(ggdata, 'vline', exact = TRUE)
  
  stopifnot(all(c('.taxa', '.rank')  %in% names(ggdata)))
  stopifnot('.taxa' %in% names(stats))
  
  if (!hasName(stats, '.adj.p'))
    stop ("No p-values to apply `p.top` constraint on.")
  
  
  
  #________________________________________________________
  # Locate the top taxa for each rank.
  #________________________________________________________
  if (!hasName(stats, '.rank'))
    stats[['.rank']] <- factor(as.character(ggdata[['.rank']][[1]]))
  
  keep_taxa <- plyr::dlply(stats, ply_cols('.rank'), function (x) {
    taxa_min_p <- split(x[['.adj.p']], x[['.taxa']]) %>%
      sapply(base::min, 1, na.rm = TRUE) %>%
      sort()
    if (p.top >= 1) { return (head(names(taxa_min_p), p.top))
    } else          { return (names(which(taxa_min_p <= p.top))) }
  })
  
  
  
  #________________________________________________________
  # Drop rows for taxa that didn't make the cut-off.
  #________________________________________________________
  for (obj_name in c('ggdata', 'stats', 'vline')) {
    
    df <- get(obj_name, inherits = FALSE)
    if (!is.data.frame(df)) next
    
    
    attrs <- attributes(df)
    
    if (!hasName(df, '.rank'))
      df[['.rank']] <- factor(names(keep_taxa)[[1]])
    
    df %<>% plyr::ddply(ply_cols('.rank'), function (x) {
      rank <- as.character(x[['.rank']][[1]])
      x[x[['.taxa']] %in% keep_taxa[[rank]],,drop=FALSE]
    }) %>% as_rbiom_tbl()
    
    df[['.taxa']] %<>% {factor(., unique(unname(unlist(keep_taxa))))}
    
    for (i in names(attrs))
      if (is.null(attr(df, i, exact = TRUE)))
        attr(df, i) <- attrs[[i]]
    
    
    assign(obj_name, df)
  }
  
  
  params$.plot_attrs[['stats']] <- stats
  attr(ggdata, 'vline')         <- vline
  params$.ggdata                <- ggdata
  
  return (invisible(params))
}

