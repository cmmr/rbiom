

#' @rdname taxa_matrix
#' @export
taxa_table <- function (
    biom, rank = -1, taxa = 6, lineage = FALSE, 
    md = ".all", unc = "singly", other = FALSE, 
    transform = "none", ties = "random", seed = 0 ) {
  
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
          transform   = transform,
          ties    = ties, 
          seed    = seed )
        
        
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
    if      (eq(transform, 'percent')) { "Relative Abundance" }
    else if (is.null(biom$depth))  { "Unrarefied Counts"  }
    else                           { "Rarefied Counts"    }
  })
  
  if (eq(params$transform, 'rank'))
    resp_label %<>% tolower() %>% paste0("Ranked Abundance\n", "(from ", ., ")")
  
  if (!params$transform %in% c('none', 'rank', 'percent'))
    resp_label %<>% paste0("\n(", params$transform, " transformed)")
  
  
  
  tbl %<>% as_rbiom_tbl()
  attr(tbl, 'response') <- ".abundance"
  attr(tbl, 'resp_label') <- resp_label
  
  
  attr(tbl, 'cmd') <- current_cmd('taxa_table')
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}




#' Taxa abundances per sample.
#' 
#' \describe{
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
#' \describe{
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
#'     hmp50$ranks
#'     
#'     taxa_matrix(hmp50, 'Phylum')[1:4,1:6]
#'     
#'     taxa_table(hmp50, 'Phylum')

taxa_matrix <- function (
    biom, rank = -1, taxa = NULL, lineage = FALSE, 
    sparse = FALSE, unc = "singly", other = FALSE, 
    transform = "none", ties = "random", seed = 0 ) {
  
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
  validate_var_choices('transform', c("none", "rank", "log", "log1p", "sqrt", "percent"))
  validate_var_choices('ties', c("average", "first", "last", "random", "max", "min"))
  validate_var_range('seed', n = 1, int = TRUE)
  
  
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
      if (eq(transform, 'percent')) counts <- rescale_cols(counts)
      mtx <- rollup(counts[names(map),], 1L, map, sum)
      mtx <- mtx[order(tolower(rownames(mtx))), colnames(counts), drop=FALSE]
      
      stopifnot(inherits(mtx, "simple_triplet_matrix"))
      
      return (mtx)
    }))
  
  
  
  #________________________________________________________
  # Only retain taxa of interest (by abundance or name).
  #________________________________________________________
  if (!is_null(taxa)) {
    
    if (is.numeric(taxa) && length(taxa) == 1) {
      rel <- sort(slam::row_means(t(t(mtx) / col_sums(mtx))), decreasing = TRUE)
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
        col_sums() %>%
        matrix(., nrow = 1, dimnames = list(other, names(.))) %>%
        as.simple_triplet_matrix() %>%
        rbind(mtx[taxa,,drop=FALSE], .)
      
      taxa <- c(taxa, other)
    }
    
  }
  
  
  
  #________________________________________________________
  # Optionally transform the computed abundance values.
  #________________________________________________________
  if (eq(transform, "rank")) {
    
    mtx %<>% as.matrix() # need those zeroes
    set.seed(seed)
    mtx[] <- base::rank(mtx, ties.method = ties)
    mtx %<>% as.simple_triplet_matrix()
    
  } else if (transform %in% c("log", "log1p", "sqrt")) {
    mtx$v <- do.call(`::`, list('base', transform))(mtx$v)
  }
  
  
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
#' @param sort  Sort the result. Options: `NULL`, `"asc"`, or `"desc"`, where
#'        `NULL` will not sort the result. `"asc"` will sort in ascending order 
#'        (smallest to largest), and `"desc"` in descending order (largest to 
#'        smallest). Ignored when the result is not a simple numeric vector. 
#'        Default: `NULL`
#' 
#' @param FUN  The function to apply to each row of the `taxa_matrix()`.
#' 
#' @param ...  Optional arguments to `FUN`.
#'        
#' @return For `taxa_sums` and `taxa_means`, a named numeric vector.
#'         For `taxa_apply`, a named vector or list with the results of `FUN`.
#'         The names are the taxa IDs.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     taxa_sums(hmp50) %>% head(4)
#'     
#'     taxa_means(hmp50, 'Family') %>% head(5)
#'     
#'     taxa_apply(hmp50, max) %>% head(5)
#'     
#'     taxa_apply(hmp50, fivenum) %>% head(5)

taxa_sums <- function (biom, rank = -1, sort = NULL, lineage = FALSE, unc = "singly", transform = "none") {
  taxa_apply(biom, FUN = base::sum, rank, sort, lineage, unc, transform)
}



#' @rdname taxa_sums
#' @export

taxa_means <- function (biom, rank = -1, sort = NULL, lineage = FALSE, unc = "singly", transform = "none") {
  taxa_apply(biom, FUN = base::mean, rank, sort, lineage, unc, transform)
}



#' @rdname taxa_sums
#' @export

taxa_apply <- function (biom, FUN, rank = -1, sort = NULL, lineage = FALSE, unc = "singly", transform = "none", ...) {
  
  stopifnot(is.function(FUN))
  validate_var_choices('sort', choices = c('asc', 'desc'), null_ok = TRUE)
  
  mtx <- taxa_matrix(biom = biom, rank = rank, lineage = lineage, sparse = TRUE, unc = unc, transform = transform)
  
  res <- if (identical(FUN, base::sum))  { slam::row_sums(mtx)
  } else if (identical(FUN, base::mean)) { slam::row_means(mtx)
  } else { slam::rowapply_simple_triplet_matrix(mtx, FUN, ...) }
  
  if (is.numeric(res) && !is.null(sort))
    res <- base::sort(res, decreasing = (sort == 'desc'))
  
  return (res)
}





#' @rdname subset
#' @export

subset_taxa <- function (x, subset, clone = TRUE, ...) {
  biom <- if (isTRUE(clone)) x$clone() else x
  df   <- biom$taxonomy
  keep <- eval(expr = substitute(subset), envir = df, enclos = parent.frame())
  df   <- biom$counts[keep,,drop = FALSE]
  suppressWarnings(biom$counts <- df)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
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
  
  
  if (!is.finite(p.top))          return (invisible(params))
  if (!hasName(stats,  '.adj.p')) return (invisible(params))
  if (!hasName(ggdata, '.taxa'))  return (invisible(params))
  if (!hasName(ggdata, '.rank'))  return (invisible(params))
  
  
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
