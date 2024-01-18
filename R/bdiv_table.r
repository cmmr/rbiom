#' Distance / dissimilarity between samples.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_default
#' 
#' @family beta_diversity
#'     
#' @param md  Include metadata in the output data frame? Options are: 
#'        \itemize{
#'          \item{`NULL` - }{ Don't include metadata. }
#'          \item{`".all"` - }{ Include all metadata. }
#'          \item{\emph{character vector} - }{
#'            Include only the specified metadata columns.
#'            Prefix the column name(s) with `==` or `!=` to limit comparisons 
#'            to within or between groups, respectively. }
#'        }
#'        Default: `NULL`
#' 
#' @param weighted  Take relative abundances into account. When 
#'        `weighted=FALSE`, only presence/absence is considered.
#'        `bdiv_table()` can accept multiple values.
#'        Default: `TRUE`
#' 
#' @param bdiv  Beta diversity distance algorithm to use. Options are:
#'        `"Bray-Curtis"`, `"Manhattan"`, `"Euclidean"`, 
#'        `"Jaccard"`, and `"UniFrac"`. `bdiv_table()` can accept multiple 
#'        values. Default: `"Bray-Curtis"`
#' 
#' @return
#' \itemize{
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
    md = NULL, within = NULL, between = NULL, trans = "none" ) {
  
  biom <- as_rbiom(biom)
  
  
  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  validate_bdiv(max = Inf)
  validate_bool('weighted', max = Inf)
  validate_meta('md', null_ok = TRUE, max = Inf, cmp = TRUE)
  validate_meta_cmp('md') # Validates and appends to `within` and `between`.
  
  
  
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
          biom     = biom, 
          bdiv     = b, 
          weighted = w, 
          tree     = tree, 
          within   = within, 
          between  = between, 
          trans    = trans )
        
        
        #________________________________________________________
        # Convert to long form
        #________________________________________________________
        tibble(
          '.sample1'  = rownames(mtx)[row(mtx)],
          '.sample2'  = colnames(mtx)[col(mtx)],
          '.weighted' = w,
          '.bdiv'     = b,
          '.distance' = as.numeric(mtx) ) %>%
          dplyr::filter(.sample1 < .sample2) %>%
          dplyr::filter(!is.na(.distance))
        
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
    if (is.numeric(map)) {
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
  
  
  lvls <- biom$metadata[['.sample']]
  tbl[['.sample1']] %<>% {factor(., levels = intersect(lvls, .))}
  tbl[['.sample2']] %<>% {factor(., levels = intersect(lvls, .))}
  
  
  attr(tbl, 'response') <- ".distance"
  
  
  return (tbl)
}





#' @rdname bdiv_table
#' @export
bdiv_matrix <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, 
    within = NULL, between = NULL, trans = "none" ) {
  
  #________________________________________________________
  # Take care not to cache filepath to tree.
  #________________________________________________________
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment())
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_bdiv()
  validate_bool("weighted")
  validate_meta('within',  null_ok = TRUE, max = Inf)
  validate_meta('between', null_ok = TRUE, max = Inf)
  validate_var_choices('trans', c("none", "rank", "log", "log1p", "sqrt"))
  
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
    
    return (x)
  })
  
  
  
  #________________________________________________________
  # Run C++ implemented dissimilarity algorithms.
  #________________________________________________________
  if (bdiv == "UniFrac") {
    
    tree      <- biom$tree
    distances <- par_unifrac(counts, pairs - 1L, tree, ifelse(weighted, 1L, 0L))
    
  } else {
    counts    <- t(as.matrix(counts))
    distances <- par_beta_div(counts, pairs - 1L, tolower(bdiv), ifelse(weighted, 1L, 0L))
  }
  
  
  
  #________________________________________________________
  # Optionally transform the computed distance values.
  #________________________________________________________
  if (trans != "none") # "rank", "log", "log1p", "sqrt"
    distances <- do.call(`::`, list('base', trans))(distances)
  
  
  
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
  set_cache_value(cache_file, mtx)
  return (mtx)
}




#' @rdname bdiv_table
#' @export
bdiv_distmat <- function (
    biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, 
    within = NULL, between = NULL, trans = "none" ) {
  
  as.dist(bdiv_matrix(
    biom     = biom, 
    bdiv     = bdiv, 
    weighted = weighted, 
    tree     = tree, 
    within   = within, 
    between  = between, 
    trans    = trans ))
}
