#' Distance / dissimilarity between samples.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' 
#' @param ...   Not used.
#' 
#' @return
#' \describe{
#'   \item{`bdiv_matrix()` - }{ An R matrix of samples x samples. }
#'   \item{`bdiv_distmat()` - }{ A dist-class distance matrix. }
#'   \item{`bdiv_table()` - }{
#'     A tibble data.frame with columns named .sample1, .sample2, 
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
#'     bdiv_table(biom, 'w_unifrac', md = ".all")
#'     
#'     # Only look at distances among the stool samples
#'     bdiv_table(biom, 'w_unifrac', md = c("==Body Site", "Sex"))
#'     
#'     # Or between males and females
#'     bdiv_table(biom, 'w_unifrac', md = c("Body Site", "!=Sex"))
#'     
#'     # All-vs-all matrix
#'     bdiv_matrix(biom, 'w_unifrac')
#'     
#'     # All-vs-all distance matrix
#'     dm <- bdiv_distmat(biom, 'w_unifrac')
#'     dm
#'     plot(hclust(dm))
#'

bdiv_table <- function (
    biom, bdiv = "bray", weighted = NULL, normalized = NULL, 
    tree = NULL, md = ".all", within = NULL, between = NULL, delta = '.all', 
    transform = "none", ties = "random", seed = 0, alpha = 0.5, cpus = NULL, ... ) {
  
  eval_dots('bdiv_table', ...)
  
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
  validate_bdiv(multiple = TRUE)
  
  validate_biom_field('md',      null_ok = TRUE, max = Inf)
  validate_biom_field('delta',   null_ok = TRUE, max = Inf)
  validate_biom_field('within',  null_ok = TRUE, max = Inf, col_type = "cat")
  validate_biom_field('between', null_ok = TRUE, max = Inf, col_type = "cat")
  
  
  
  
  #________________________________________________________
  # Multiple combinations of bdiv into one table.
  #________________________________________________________
  tbl <- NULL
  
  for (i in seq_along(bdiv))
    tbl %<>% dplyr::bind_rows(local({
      
      #________________________________________________________
      # Compute the distance matrix
      #________________________________________________________
      mtx <- bdiv_matrix(
        biom      = biom, 
        bdiv      = bdiv[[i]], 
        tree      = tree, 
        within    = within, 
        between   = between, 
        transform = transform, 
        ties      = ties, 
        seed      = seed, 
        alpha     = alpha, 
        cpus      = cpus )
      
      
      #________________________________________________________
      # Convert to long form
      #________________________________________________________
      lt <- as.vector(lower.tri(mtx))
      tibble(
        '.sample1'  = colnames(mtx)[col(mtx)][lt],
        '.sample2'  = rownames(mtx)[row(mtx)][lt],
        '.bdiv'     = bdiv[[i]],
        '.distance' = as.vector(mtx)[lt] ) %>%
        dplyr::filter(!is.na(.data$.distance))
      
    }))
  
  tbl[['.bdiv']] %<>% factor(levels = unique(bdiv))
  
  
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
  
  resp_label <- "Beta Dissimilarity"
  
  if (length(bdiv) == 1)
    resp_label <- ecodive::match_metric(bdiv)$name
  
  if (eq(params$transform, 'rank'))
    resp_label <- paste("Ranked", resp_label)
  
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
    biom, bdiv = "bray", weighted = NULL, normalized = NULL, 
    tree = NULL, within = NULL, between = NULL, transform = "none", 
    ties = "random", seed = 0, alpha = 0.5, cpus = NULL ) {
  
  params <- eval_envir(environment())
  cmd    <- sprintf("bdiv_matrix(%s)", as.args(params, fun = bdiv_matrix))
  remove("params")
  
  validate_bdiv()
  
  mtx <- as.matrix(bdiv_distmat(
    biom      = biom,
    bdiv      = bdiv,
    tree      = tree,
    within    = within,
    between   = between,
    transform = transform,
    ties      = ties,
    seed      = seed,
    alpha     = alpha, 
    cpus      = cpus ))
  
  attr(mtx, 'cmd') <- cmd
  return (mtx)
}




#' @rdname bdiv_table
#' @export
bdiv_distmat <- function (
    biom, bdiv = "bray", weighted = NULL, normalized = NULL, alpha = 0.5, 
    tree = NULL, within = NULL, between = NULL, transform = "none", 
    ties = "random", seed = 0, cpus = NULL ) {
  
  #________________________________________________________
  # Take care not to cache filepath to tree.
  #________________________________________________________
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment())
  cmd    <- sprintf("bdiv_distmat(%s)", as.args(params, fun = bdiv_matrix))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_distmat', params)
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
  validate_alpha()
  validate_var_choices('transform', c("none", "rank", "log", "log1p", "sqrt", "percent"))
  validate_var_choices('ties', c("average", "first", "last", "random", "max", "min"))
  validate_seed()
  validate_cpus()
  
  
  
  #________________________________________________________
  # Only calculate distances between relevant samples.
  #________________________________________________________
  pairs <- local({
    
    if (length(c(within, between)) == 0)
      return (NULL)
    
    # All unique pairing combinations.
    sids <- biom$samples
    x    <- combn(factor(sids, levels = sids), 2)
    i    <- replicate(ncol(x), TRUE)
    
    # Limit to only within or between comparisons.
    for (col in c(within, between)) {
      map <- pull(biom, col)
      if (col %in% within)  i <- i & (map[x[1,]] == map[x[2,]])
      if (col %in% between) i <- i & (map[x[1,]] != map[x[2,]])
    }
    
    return (which(i))
  })
  
  
  
  #________________________________________________________
  # Run dissimilarity algorithms implemented in C.
  #________________________________________________________
  dm <- ecodive::beta_div(
    counts = biom, 
    metric = bdiv, 
    alpha  = alpha, 
    tree   = if (is.null(tree)) biom$tree else tree,
    pairs  = pairs, 
    cpus   = cpus )
  
  
  #________________________________________________________
  # Optionally transform the computed distance values.
  #________________________________________________________
  if (eq(transform, "rank")) {
    
    if (eq(ties, 'random')) { # Preserve current .Random.seed
      oldseed <- if (exists(".Random.seed")) .Random.seed else NULL
      set.seed(seed)
      if (!is.null(oldseed)) on.exit(.Random.seed <- oldseed)
    }
    dm[] <- base::rank(dm, na.last = 'keep', ties.method = ties)
    
  } else if (eq(transform, "percent")) {
    dm <- tryCatch(dm / sum(dm), error = function (e) { warning(e); dm })
    
  } else if (transform %in% c("log", "log1p", "sqrt")) {
    dm <- do.call(`::`, list('base', transform))(dm)
  }
  
  
  
  #________________________________________________________
  # Return the distance matrix
  #________________________________________________________
  
  attr(dm, 'cmd')    <- cmd
  attr(dm, 'metric') <- bdiv
  set_cache_value(cache_file, dm)
}
