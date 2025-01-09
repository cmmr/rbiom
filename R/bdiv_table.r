#' Distance / dissimilarity between samples.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' 
#' @return
#' \describe{
#'   \item{`bdiv_matrix()` - }{ An R matrix of samples x samples. }
#'   \item{`bdiv_distmat()` - }{ A dist-class distance matrix. }
#'   \item{`bdiv_table()` - }{
#'     A tibble data.frame with columns names .sample1, .sample2, .weighted, 
#'     .bdiv, .distance, and any fields requested by `md`. Numeric metadata 
#'     fields will be returned as `abs(x - y)`; categorical metadata fields as 
#'     `"x"`, `"y"`, or `"x vs y"`. }
#' }
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Subset to four samples
#'     biom <- hmp50$clone()
#'     biom$counts <- biom$counts[,c("HMP18", "HMP19", "HMP20", "HMP21")]
#'     
#'     # Return in long format with metadata
#'     bdiv_table(biom, 'unifrac', md = ".all")
#'     
#'     # Only look at distances among the stool samples
#'     bdiv_table(biom, 'unifrac', md = c("==Body Site", "Sex"))
#'     
#'     # Or between males and females
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "!=Sex"))
#'     
#'     # All-vs-all matrix
#'     bdiv_matrix(biom, 'unifrac')
#'     
#'     # All-vs-all distance matrix
#'     dm <- bdiv_distmat(biom, 'unifrac')
#'     dm
#'     plot(hclust(dm))
#'

bdiv_table <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, 
    md = ".all", within = NULL, between = NULL, delta = '.all', 
    transform = "none", ties = "random", seed = 0, cpus = -1 ) {
  
  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- sprintf("bdiv_table(%s)", as.args(params, fun = bdiv_table))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_table', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Strips '==' and '!='; appends to within and between.
  #________________________________________________________
  validate_var_cmp(c('md', 'delta'))
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  validate_bdiv(max = Inf)
  validate_bool('weighted', max = Inf)
  
  validate_biom_field('md',      null_ok = TRUE, max = Inf)
  validate_biom_field('delta',   null_ok = TRUE, max = Inf)
  validate_biom_field('within',  null_ok = TRUE, max = Inf, col_type = "cat")
  validate_biom_field('between', null_ok = TRUE, max = Inf, col_type = "cat")
  
  
  
  
  #________________________________________________________
  # Multiple combinations of bdiv/weighted into one table.
  #________________________________________________________
  tbl <- NULL
  
  for (w in weighted)
    for (b in bdiv)
      tbl %<>% dplyr::bind_rows(local({
        
        #________________________________________________________
        # Compute the distance matrix
        #________________________________________________________
        mtx <- bdiv_matrix(
          biom      = biom, 
          bdiv      = b, 
          weighted  = w, 
          tree      = tree, 
          within    = within, 
          between   = between, 
          transform = transform, 
          ties      = ties, 
          seed      = seed, 
          cpus      = cpus )
        
        
        #________________________________________________________
        # Convert to long form
        #________________________________________________________
        tibble(
          '.sample1'  = rownames(mtx)[row(mtx)],
          '.sample2'  = colnames(mtx)[col(mtx)],
          '.weighted' = w,
          '.bdiv'     = b,
          '.distance' = as.numeric(mtx) ) %>%
          dplyr::filter(.data$.sample1 < .data$.sample2) %>%
          dplyr::filter(!is.na(.data$.distance))
        
      }))
  
  tbl[['.bdiv']] %<>% factor(levels = bdiv)
  
  
  #________________________________________________________
  # Add metadata columns
  #________________________________________________________
  for (col in unique(c(md, within, between))) {
    
    map <- pull(biom, col)
    v1  <- map[tbl[['.sample1']]]
    v2  <- map[tbl[['.sample2']]]
    
    
    # Compute abs(v1 - v2)
    #________________________________________________________
    if (is.numeric(map) && col %in% delta) {
      tbl[[col]] <- abs(v1 - v2)
    }
    
    
    # Transform to "v1 vs v2"
    #________________________________________________________
    else {
      
      
      # Change "Male vs Male" --> "Male", and sort by factor
      # level, e.g. "Male vs Female" --> "Female vs Male".
      #________________________________________________________
      tbl[[col]] <- ifelse(
        test = (v1 == v2), 
        yes  = as.character(v1), 
        no   = ifelse(
          test = as.numeric(v1) < as.numeric(v2), 
          yes  = paste(v1, "vs", v2), 
          no   = paste(v2, "vs", v1) ))
      
      
      # Keep factors as factors.
      #________________________________________________________
      if (is.factor(map)) {
        
        lvls <- levels(map)
        
        if (length(lvls) > 1)
          lvls %<>% {c(., apply(combn(., 2), 2L, paste, collapse=" vs "))}
        
        tbl[[col]] %<>% {factor(., levels = intersect(lvls, .))}
      }
      
    }
    
    
  }
  
  
  
  #________________________________________________________
  # Descriptive label for y-axis.
  #________________________________________________________
  resp_label <- local({
    
    if (length(bdiv) == 1 && length(weighted) == 1) {
      w <- ifelse(weighted, "Weighted", "Unweighted")
      bdiv %<>% paste("Distance")
      
      if (eq(params$transform, 'rank')) { paste0("Ranked ", w, "\n", bdiv) }
      else                          { paste(w, bdiv)                   }
      
    } else {
      if (eq(params$transform, 'rank')) { "Ranked Beta Dissimilarity" }
      else                          { "Beta Dissimilarity"        }
    }
    
  })
  
  if (!params$transform %in% c('none', 'rank', 'percent'))
    resp_label %<>% paste0("\n(", params$transform, " transformed)")
  
  
  
  tbl %<>% as_rbiom_tbl()
  attr(tbl, 'cmd')      <- cmd
  attr(tbl, 'response') <- ".distance"
  attr(tbl, 'resp_label') <- resp_label
  
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}





#' @rdname bdiv_table
#' @export
bdiv_matrix <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, 
    within = NULL, between = NULL, 
    transform = "none", ties = "random", seed = 0, cpus = -1 ) {
  
  #________________________________________________________
  # Take care not to cache filepath to tree.
  #________________________________________________________
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment())
  cmd    <- sprintf("bdiv_matrix(%s)", as.args(params, fun = bdiv_matrix))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_matrix', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  #________________________________________________________
  # Ensure `within` and `between` are non-overlapping.
  #________________________________________________________
  validate_biom_field('within',  null_ok = TRUE, max = Inf, col_type = "cat")
  validate_biom_field('between', null_ok = TRUE, max = Inf, col_type = "cat")
  
  if (length(x <- intersect(within, between)) > 0)
    cli_abort("Metadata field name{?s} {.val {x}} cannot be set as both a within (==) and between (!=) grouping.")
  
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_bdiv()
  validate_bool("weighted")
  validate_var_choices('transform', c("none", "rank", "log", "log1p", "sqrt", "percent"))
  validate_var_choices('ties', c("average", "first", "last", "random", "max", "min"))
  validate_var_range('seed', n = 1, int = TRUE)
  validate_var_range('cpus', n = 1, int = TRUE, range = c(-1, 256))
  
  if (!is.null(tree)) {
    biom <- biom$clone()
    biom$tree <- tree
  }
  
  
  
  #________________________________________________________
  # Order the sparse matrix's values by sample, then by taxa
  #________________________________________________________
  if (bdiv == "UniFrac") {
    
    if (is.null(biom$tree))
      cli_abort(c(x = "Phylogenetic tree required for UniFrac."))
    
    counts <- biom$counts
    counts <- counts[as.character(biom$tree$tip.label),]
    
  } else {
    counts <- biom$counts
  }
  
  ord      <- order(counts$j, counts$i)
  counts$i <- counts$i[ord]
  counts$j <- counts$j[ord]
  counts$v <- counts$v[ord]
  
  
  
  
  #________________________________________________________
  # Only calculate distances between relevant samples.
  #________________________________________________________
  sids  <- colnames(counts)
  pairs <- local({
    
    # Limit to sample1 < sample2.
    x <- combn(factor(sids, levels = sids), 2)
    
    # Limit to only within or between comparisons.
    for (col in c(within, between)) {
      map <- pull(biom, col)
      if (col %in% within)  x <- x[,map[x[1,]] == map[x[2,]]]
      if (col %in% between) x <- x[,map[x[1,]] != map[x[2,]]]
    }
    
    x <- apply(x, 1L, as.integer)
    
    if (!is.matrix(x))
      x <- matrix(x, nrow = 1)
    
    return (x)
  })
  
  
  
  #________________________________________________________
  # Run C++ implemented dissimilarity algorithms.
  #________________________________________________________
  if (bdiv == "UniFrac") {
    
    tree      <- biom$tree
    distances <- parallel_unifrac(counts, pairs - 1L, tree, weighted, cpus)
    
  } else {
    
    bdiv_func  <- switch(
      EXPR = bdiv,
      'Jaccard'     = if (weighted) bdiv_jaccard_w    else bdiv_jaccard_u, 
      'Bray-Curtis' = if (weighted) bdiv_braycurtis_w else bdiv_braycurtis_u, 
      'Manhattan'   = if (weighted) bdiv_manhattan_w  else bdiv_manhattan_u, 
      'Euclidean'   = if (weighted) bdiv_euclidean_w  else bdiv_euclidean_u )
    
    vec_counts <- as.vector(as.matrix(counts))
    vec_pairs  <- as.vector(t(pairs) - 1L)
    distances  <- bdiv_func(vec_counts, nrow(counts), vec_pairs, nrow(pairs), cpus)
  }
  
  
  
  #________________________________________________________
  # Optionally transform the computed distance values.
  #________________________________________________________
  if (eq(transform, "rank")) {
    set.seed(seed)
    distances <- base::rank(distances, ties.method = ties)
    
  } else if (eq(transform, "percent")) {
    distances <- tryCatch(distances / sum(distances), error = function (e) { warning(e); distances })
    
  } else if (transform %in% c("log", "log1p", "sqrt")) {
    distances <- do.call(`::`, list('base', transform))(distances)
  }
  
  
  
  #________________________________________________________
  # Assemble a complete all-vs-all matrix.
  #________________________________________________________
  n   <- length(sids)
  mtx <- matrix(data = rep(as.numeric(NA), n * n), nrow = n)
  mtx[pairs]          <- distances
  mtx[pairs[,c(2,1)]] <- distances
  mtx[cbind(1:n,1:n)] <- rep(0, n)
  dimnames(mtx)       <- list(sids, sids)
  
  
  
  #________________________________________________________
  # Return the matrix
  #________________________________________________________
  attr(mtx, 'cmd') <- cmd
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}




#' @rdname bdiv_table
#' @export
bdiv_distmat <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, 
    within = NULL, between = NULL, transform = "none", cpus = -1 ) {
  
  stats::as.dist(bdiv_matrix(
    biom      = biom, 
    bdiv      = bdiv, 
    weighted  = weighted, 
    tree      = tree, 
    within    = within, 
    between   = between, 
    transform = transform,
    cpus      = cpus ))
}
