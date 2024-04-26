

#' @rdname taxa_matrix
#' @export
taxa_table <- function (
    biom, rank = -1, taxa = 6, lineage = FALSE, 
    md = ".all", unc = "singly", other = FALSE, trans = "none" ) {
  
  biom <- as_rbiom(biom)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env()
  cache_file <- get_cache_file('taxa_table', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  params <- list2env(params)
  tbl <- with(params, {
    
    
    #________________________________________________________
    # Validate user's arguments.
    #________________________________________________________
    validate_rank(max = Inf)
    validate_biom_field('md', max = Inf, null_ok = TRUE)
    
    
    
    #________________________________________________________
    # Return multiple ranks in a single table.
    #________________________________________________________
    tbl <- NULL
    
    for (r in rank)
      tbl %<>% dplyr::bind_rows(local({
        
        mtx <- taxa_matrix(
          biom    = biom, 
          rank    = r, 
          taxa    = taxa, 
          lineage = lineage, 
          sparse  = FALSE, 
          unc     = unc, 
          other   = other,
          trans   = trans )
        
        
        #________________________________________________________
        # Pivot Longer
        #________________________________________________________
        tibble(
          '.rank'      = r,
          '.sample'    = colnames(mtx)[col(mtx)],
          '.taxa'      = rownames(mtx)[row(mtx)],
          '.abundance' = as.numeric(mtx) )
        
      }))
    
    tbl[['.rank']]   %<>%  factor(., levels = rank)
    tbl[['.sample']] %<>% {factor(., levels = intersect(biom$samples, .))}
    tbl[['.taxa']]   %<>% {factor(., levels = unique(.))}
    
    
    
    #________________________________________________________
    # Add Metadata
    #________________________________________________________
    if (length(md) > 0)
      tbl %<>% left_join( 
        by = '.sample',
        y  = biom$metadata[,unique(c('.sample', md))] )
    
    tbl
  })
  
  
  
  #________________________________________________________
  # Descriptive label for y-axis.
  #________________________________________________________
  resp_label <- with(params, {
    if      (eq(trans, 'percent')) { "Relative Abundance" }
    else if (is.null(biom$depth))  { "Unrarefied Counts"  }
    else                           { "Rarefied Counts"    }
  })
  
  
  
  tbl %<>% as_rbiom_tbl()
  attr(tbl, 'response') <- ".abundance"
  attr(tbl, 'resp_label') <- resp_label
  
  
  attr(tbl, 'cmd') <- current_cmd('taxa_table')
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}




#' Taxa abundances per sample.
#' 
#' \itemize{
#'   \item{`taxa_matrix()` - }{ Accepts a single `rank` and returns a matrix. }
#'   \item{`taxa_table()` - }{ Can accept more than one `rank` and returns a tibble data.frame.  }
#' }
#' 
#' 
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' 
#' @return 
#' \itemize{
#'   \item{`taxa_matrix()` - }{
#'     A numeric matrix with taxa as rows, and samples as columns. }
#'   \item{`taxa_table()` - }{
#'     A tibble data frame with column names .sample, .taxa, .abundance, and any requested by `md`. }
#' }
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom$ranks
#'     
#'     taxa_matrix(hmp50, 'Phylum')[1:4,1:6]
#'     
#'     taxa_table(hmp50, 'Phylum')

taxa_matrix <- function (
    biom, rank = -1, taxa = NULL, lineage = FALSE, 
    sparse = FALSE, unc = "singly", other = FALSE, trans = "none" ) {
  
  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- sprintf("taxa_matrix(%s)", as.args(params, fun = taxa_matrix))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_matrix', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_taxa(null_ok = TRUE)
  validate_rank()
  validate_bool("sparse")
  validate_var_choices('trans', c("none", "rank", "log", "log1p", "sqrt", "percent"))
  
  
  #________________________________________________________
  # Fetch taxonomy map with ambiguous names corrected.
  #________________________________________________________
  map <- taxa_map(
    biom    = biom, 
    rank    = rank, 
    unc     = unc, 
    lineage = lineage )
  
  
  
  #________________________________________________________
  # Compute abundance matrix for this rank.
  #________________________________________________________
  mtx <- tryCatch(
    error = function (e) stop("Unable to group by taxonomic level: ", e),
    expr  = local({
      
      counts <- biom$counts
      if (eq(trans, 'percent')) counts <- rescale_cols(counts)
      mtx <- slam::rollup(counts[names(map),], 1L, map, sum)
      mtx <- mtx[order(tolower(rownames(mtx))), colnames(counts), drop=FALSE]
      
      stopifnot(is(mtx, "simple_triplet_matrix"))
      
      return (mtx)
    }))
  
  
  
  #________________________________________________________
  # Only retain taxa of interest (by abundance or name).
  #________________________________________________________
  if (!is_null(taxa)) {
    
    if (is.numeric(taxa) && length(taxa) == 1) {
      rel <- sort(slam::row_means(t(t(mtx) / slam::col_sums(mtx))), decreasing = TRUE)
      if (taxa >= 1) { taxa <- head(names(rel), taxa) 
      } else         { taxa <- names(rel)[rel >= taxa] }
      
    } else if (is.character(taxa)) {
      taxa <- intersect(as.character(taxa), rownames(mtx))
      
    } else {
      stop ("Invalid argument for `taxa`: ", capture.output(str(taxa)))
    }
    
    
    if (length(taxa) == 0)
      stop("No taxa match the criteria: ", capture.output(str(taxa)))
    
    
    if (isFALSE(other) || is.null(other)) {
      
      mtx <- mtx[taxa,,drop=FALSE]
      
    } else {
      
      if (!is_scalar_character(other) || is_na(other))
        other <- "Other"
      
      mtx <- mtx[setdiff(rownames(mtx), taxa),,drop=FALSE] %>% 
        slam::col_sums() %>%
        matrix(., nrow = 1, dimnames = list(other, names(.))) %>%
        slam::as.simple_triplet_matrix() %>%
        rbind(mtx[taxa,,drop=FALSE], .)
      
      taxa <- c(taxa, other)
    }
    
  }
  
  
  
  #________________________________________________________
  # Optionally transform the computed abundance values.
  #________________________________________________________
  if (trans %in% c("rank", "log", "log1p", "sqrt"))
    mtx$v <- do.call(`::`, list('base', trans))(mtx$v)
  
  
  #________________________________________________________
  # Convert from slam to base matrix if sparse=FALSE
  #________________________________________________________
  if (isFALSE(sparse))
    mtx <- as.matrix(mtx)
  
  
  attr(mtx, 'cmd') <- cmd
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}





#' Get summary taxa abundances.
#' 
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' 
#' @param rank  The taxonomic rank to return sums or means for. The default, 
#'        \code{0}, returns per-OTU summaries.
#'        
#' @return A named, sorted numeric vector.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     taxa_sums(hmp50) %>% head(4)
#'     
#'     taxa_means(hmp50, 'Family') %>% head(5)
#'

taxa_sums <- function (biom, rank = 0) {
  
  taxa_matrix(biom, rank, sparse = TRUE) %>%
    slam::row_sums() %>%
    sort(decreasing = TRUE)
}



#' @rdname taxa_sums
#' @export

taxa_means <- function (biom, rank = 0) {
  
  taxa_matrix(biom, rank, sparse = TRUE) %>%
    slam::row_means() %>%
    sort(decreasing = TRUE)
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
