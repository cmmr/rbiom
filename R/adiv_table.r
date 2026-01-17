
#' Calculate the alpha diversity of each sample.
#' 
#' @inherit documentation_default
#' 
#' @family alpha_diversity
#'        
#' @return 
#' \describe{
#'   \item{`adiv_vector()` - }{ A named numeric vector. }
#'   \item{`adiv_matrix()` - }{
#'     A matrix of samples x metric. 
#'     The first column, 'depth', is never transformed. }
#'   \item{`adiv_table()` - }{
#'     A tibble data.frame of alpha diversity values.
#'     Each combination of sample/\code{adiv} has its own row.
#'     Column names are \bold{.sample}, \bold{.depth}, \bold{.adiv}, 
#'     and \bold{.diversity}, followed by any metadata fields requested by 
#'     \code{md}. }
#' }
#' 
#' @seealso [sample_sums()] for sample depths.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- hmp50[1:5]
#'     
#'     adiv_table(biom)
#'     
#'     biom <- rarefy(biom)
#'     adiv_table(biom, md = NULL)
#'     
#'     adiv_vector(biom, 'faith')
#'     
#'     adiv_matrix(biom)

adiv_table <- function (
    biom, adiv = "shannon", md = ".all", tree = NULL, transform = "none", 
    ties = "random", seed = 0, cpus = NULL ) {
  
  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- sprintf("adiv_table(%s)", as.args(params, fun = adiv_table))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('adiv_table', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_adiv(multiple = TRUE)
  validate_biom_field('md', null_ok = TRUE, max = Inf)
  
  
  #________________________________________________________
  # Compute adiv values
  #________________________________________________________
  mtx <- adiv_matrix(
    biom      = biom, 
    adiv      = adiv, 
    tree      = tree, 
    transform = transform,
    ties      = ties,
    seed      = seed,
    cpus      = cpus )
  
  tbl <- tibble(
    '.sample'    = rownames(mtx)[row(mtx)] %>% factor(levels = rownames(mtx)),
    '.depth'     = unname(mtx[,'depth'][row(mtx)]),
    '.adiv'      = colnames(mtx)[col(mtx)] %>% factor(levels = adiv),
    '.diversity' = as.numeric(mtx) ) %>%
    dplyr::filter(!is.na(.data$.adiv))
  
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (!is.null(md))
    tbl %<>% left_join(
      y  = biom$metadata[,unique(c('.sample', md))], 
      by = '.sample' )
  
  
  
  #________________________________________________________
  # Descriptive label for y-axis.
  #________________________________________________________
  if (length(adiv) == 1) { resp_label <- ecodive::match_metric(adiv)$name } 
  else                   { resp_label <- "Alpha Diversity"                }
  
  
  
  tbl %<>% as_rbiom_tbl()
  attr(tbl, 'cmd')        <- cmd
  attr(tbl, 'response')   <- ".diversity"
  attr(tbl, 'resp_label') <- resp_label
  
  
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}



#' @rdname adiv_table
#' @export
adiv_matrix <- function (
    biom, adiv = c('observed', 'shannon', 'simpson'), tree = NULL, 
    transform = "none", ties = "random", seed = 0, cpus = NULL ) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment())
  cmd    <- sprintf("adiv_matrix(%s)", as.args(params, fun = adiv_matrix))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('adiv_matrix', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_adiv(multiple = TRUE)
  
  
  
  #________________________________________________________
  # We want a numeric matrix of samples x adiv metrics
  #________________________________________________________
  
  mtx <- matrix(
    data     = NA_real_, 
    nrow     = biom$n_samples, 
    ncol     = length(adiv) + 1,
    dimnames = list(biom$samples, c('depth', adiv)) )
  
  mtx[,'depth'] <- colSums(biom$counts)
  
  for (i in seq_along(adiv))
    mtx[,i+1] <- adiv_vector(
      biom      = biom, 
      adiv      = adiv[[i]], 
      tree      = tree, 
      transform = transform,
      ties      = ties,
      seed      = seed,
      cpus      = cpus )
  
  
  #________________________________________________________
  # Return the matrix of diversity values
  #________________________________________________________
  attr(mtx, 'cmd') <- cmd
  set_cache_value(cache_file, mtx)
  return (mtx)
}



#' @rdname adiv_table
#' @export
adiv_vector <- function (
    biom, adiv = "shannon", tree = NULL, transform = "none", 
    ties = "random", seed = 0, cpus = NULL ) {
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment())
  cmd    <- sprintf("adiv_vector(%s)", as.args(params, fun = adiv_vector))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('adiv_vector', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_adiv()
  validate_var_choices('transform', c("none", "rank", "log", "log1p", "sqrt"))
  validate_var_choices('ties', c("average", "first", "last", "random", "max", "min"))
  validate_seed()
  validate_cpus()
  
  
  #________________________________________________________
  # Run alpha diversity algorithms implemented in C.
  #________________________________________________________
  vec <- ecodive::alpha_div(
    counts = biom, 
    metric = adiv, 
    tree   = if (is.null(tree)) biom$tree else tree,
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
    vec <- base::rank(vec, na.last = 'keep', ties.method = ties)
    
  } else if (eq(transform, "percent")) {
    vec <- tryCatch(vec / sum(vec), error = function (e) { warning(e); vec })
    
  } else if (transform %in% c("log", "log1p", "sqrt")) {
    vec <- do.call(`::`, list('base', transform))(vec)
  }
  
  
  
  #________________________________________________________
  # Return the diversity values
  #________________________________________________________
  attr(vec, 'cmd')    <- cmd
  attr(vec, 'metric') <- adiv
  set_cache_value(cache_file, vec)
}



#' Summarize the taxa observations in each sample.
#' 
#' @inherit documentation_default
#' 
#' @family samples
#' @family rarefaction
#' @family taxa_abundance
#' 
#' @param sort  Sort the result. Options: `NULL` - don't sort; `"asc"` - in 
#'        ascending order (smallest to largest); `"desc"` - in descending order 
#'        (largest to smallest). Ignored when the result is not a simple 
#'        numeric vector. Default: `NULL`
#' 
#' @param FUN  The function to apply to each column of `taxa_matrix()`.
#' 
#' @param ...  Optional arguments to `FUN`.
#' 
#' @return For `sample_sums`, A named numeric vector of the number of 
#'         observations in each sample. For `sample_apply`, a named vector or 
#'         list with the results of `FUN`. The names are the taxa IDs.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     library(ggplot2)
#'     
#'     sample_sums(hmp50, sort = 'asc') %>% head()
#'     
#'     # Unique OTUs and "cultured" classes per sample
#'     nnz <- function (x) sum(x > 0) # number of non-zeroes
#'     sample_apply(hmp50, nnz, 'otu') %>% head()
#'     sample_apply(hmp50, nnz, 'class', unc = 'drop') %>% head()
#'     
#'     # Number of reads in each sample's most abundant family
#'     sample_apply(hmp50, base::max, 'f', sort = 'desc') %>% head()
#'     
#'     ggplot() + geom_histogram(aes(x=sample_sums(hmp50)), bins = 20)

sample_sums <- function (biom, rank = -1, sort = NULL, unc = "singly") {
  sample_apply(biom, FUN = base::sum, rank, sort, unc)
}



#' @rdname sample_sums
#' @export

sample_apply <- function (biom, FUN, rank = -1, sort = NULL, unc = "singly", ...) {
  
  stopifnot(is.function(FUN))
  validate_var_choices('sort', choices = c('asc', 'desc'), null_ok = TRUE)
  
  mtx <- taxa_matrix(biom = biom, rank = rank, sparse = TRUE, unc = unc)
  
  res <- if (identical(FUN, base::sum))  { colSums(mtx)
  } else if (identical(FUN, base::mean)) { colMeans(mtx)
  } else { apply(mtx, 2, FUN, ...) }
  
  if (is.numeric(res) && !is.null(sort))
    res <- base::sort(res, decreasing = (sort == 'desc'))
  
  return (res)
}
