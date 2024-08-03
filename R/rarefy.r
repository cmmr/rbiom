
#' Rarefy OTU counts.
#' 
#' Sub-sample OTU observations such that all samples have an equal number.
#' If called on data with non-integer abundances, values will be re-scaled to 
#' integers between 1 and `depth` such that they sum to `depth`.
#'
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' 
#' @family rarefaction
#' @family transformations
#' 
#' @param depth   How many observations to keep per sample. When 
#'        `0 < depth < 1`, it is taken as the minimum percentage of the 
#'        dataset's observations to keep. Ignored when `n` is specified.
#'        Default: `0.1`
#'
#' @param n   The number of samples to keep. When `0 < n < 1`, it is taken as 
#'        the percentage of samples to keep. If negative, that number or 
#'        percentage of samples is dropped. If `0`, all samples are kept. If 
#'        `NULL`, `depth` is used instead.
#'        Default: `NULL`
#'     
#' @param seed   An integer seed for randomizing which observations to keep or 
#'        drop. If you need to create different random rarefactions of the same 
#'        data, set the seed to a different number each time.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     sample_sums(hmp50) %>% head()
#'     
#'     biom <- rarefy(hmp50)
#'     sample_sums(biom) %>% head()
#' 
rarefy <- function (biom, depth = 0.1, n = NULL, seed = 0, clone = TRUE) {
  biom <- as_rbiom(biom)
  if (isTRUE(clone)) biom <- biom$clone()
  biom$counts <- rarefy_cols(mtx = biom$counts, depth = depth, n = n, seed = seed)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}




#' Transform a counts matrix.
#' 
#' Rarefaction subset counts so that all samples have the same number of 
#' observations. Rescaling rows or cols scales the matrix values so that row 
#' sums or column sums equal 1.
#'
#' @inherit documentation_default
#' 
#' @family rarefaction
#' @family transformations
#'     
#' @param depth   How many observations to keep per sample. When 
#'        `0 < depth < 1`, it is taken as the minimum percentage of the 
#'        dataset's observations to keep. Ignored when `n` is specified.
#'        Default: `0.1`
#'
#' @param n   The number of samples to keep. When `0 < n < 1`, it is taken as 
#'        the percentage of samples to keep. If negative, that number or 
#'        percentage of samples is dropped. If `0`, all samples are kept. If 
#'        `NULL`, `depth` is used instead.
#'        Default: `NULL`
#'     
#' @param seed   An integer to use for seeding the random number generator. If
#'        you need to create different random rarefactions of the same 
#'        matrix, set this seed value to a different number each time.
#'     
#' @return The rarefied or rescaled matrix.
#' 
#' @examples
#'     library(rbiom)
#'     
#'     # rarefy_cols --------------------------------------
#'     biom <- hmp50$clone()
#'     sample_sums(biom) %>% head(10)
#'
#'     biom$counts %<>% rarefy_cols(depth=1000)
#'     sample_sums(biom) %>% head(10)
#'     
#'     
#'     # rescaling ----------------------------------------
#'     mtx <- matrix(sample(1:20), nrow=4)
#'     mtx
#'     
#'     rowSums(mtx)
#'     rowSums(rescale_rows(mtx))
#'     
#'     colSums(mtx)
#'     colSums(rescale_cols(mtx))
#'
rarefy_cols <- function (mtx, depth = 0.1, n = NULL, seed = 0) {
  
  params <- eval_envir(environment())
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('rarefy_cols', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity checks.
  #________________________________________________________
  mtx <- as.simple_triplet_matrix(mtx)
  
  if (!(is.numeric(depth) && length(depth) == 1 && !is.na(depth)))
    cli_abort(c('x' = "{.var depth} must be a single number, not {.type {depth}}."))
  
  if (!((depth > 0 && depth < 1) || (depth >= 1 && depth %% 1 == 0)))
    cli_abort(c('x' = "{.var depth} must between 0 and 1, or a positive integer, not {depth}."))
  
  if (!is.null(n)) {
    if (!(is.numeric(n) && length(n) == 1 && !is.na(n)))
      cli_abort(c('x' = "{.var n} must be a single number, not {.type {n}}."))
    
    if (abs(n) >= 1 && n %% 1 != 0)
      cli_abort(c('x' = "{.var n} must between -1 and 1, or an integer, not {n}."))
  }
  
  if (!is_scalar_integerish(seed) || is_na(seed))
    cli_abort(c('x' = "{.var seed} must be a single integer, not {.type {seed}}."))
  
  
  
  #________________________________________________________
  # Integer data. Randomly select observations to keep.
  #________________________________________________________
  if (all(mtx$v %% 1 == 0)) {
    
    
    # Set depth according to number/pct of samples to keep/drop.
    if (!is.null(n)) {
      if (n == 0)     n <- mtx$ncol     # Keep all
      if (abs(n) < 1) n <- n * mtx$ncol # Keep/drop percentage
      if (n < -1)     n <- mtx$ncol + n # Drop n
      n     <- max(1, floor(n))         # Keep at least one
      depth <- unname(head(tail(sort(col_sums(mtx)), n), 1))
    }
    
    
    # Depth is given as minimum percent of obs. to keep.
    if (depth < 1) {
      sums  <- col_sums(mtx)
      depth <- (sum(sums) * depth) / length(sums)
      depth <- min(sums[sums >= depth])
    }
    
    mtx <- rcpp_rarefy(mtx, depth, seed)
    
  }
  
  
  #________________________________________________________
  # Rescale fractional data to integers. Total = depth.
  #________________________________________________________
  else {
    
    if (!is_scalar_integerish(depth))
      depth <- 10000
    
    mtx <- as.matrix(mtx)
    
    # Does this, but then nudges up and down until == depth:
    # round(t((t(mtx) / colSums(mtx))) * depth)
    
    mtx <- apply(mtx, 2L, function (x) { 
      x        <- x / sum(x)
      y        <- round(x * depth)
      maxIters <- length(x)
      
      while (sum(y) < depth && maxIters > 0) {
        i        <- which.min( ((y + 1) / depth) - x )
        y[i]     <- y[i] + 1
        maxIters <- maxIters - 1
      }
      while (sum(y) > depth && maxIters > 0) {
        i        <- which.min( x - ((y - 1) / depth) )
        y[i]     <- y[i] - 1
        maxIters <- maxIters - 1
      }
      
      return (y)
    })
    
  }
  
  
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}




#' @rdname rarefy_cols
#' @export
rescale_cols <- function (mtx) {
  
  if (is.simple_triplet_matrix(mtx))
    return (t(t(mtx) / col_sums(mtx)))
  
  mtx <- as.matrix(mtx)
  t(t(mtx) / colSums(mtx))
}



#' @rdname rarefy_cols
#' @export
rescale_rows <- function (mtx) {
  
  if (is.simple_triplet_matrix(mtx))
    return (mtx / col_sums(t(mtx)))
  
  mtx <- as.matrix(mtx)
  mtx / colSums(t(mtx))
}




#' Suggest a 'good' rarefaction depth.
#' 
#' @noRd
#' @keywords internal
#'     
#' @return An integer.
#'

rare_suggest <- function (mtx) {
  
  stopifnot(is.simple_triplet_matrix(mtx))
  
  if (all(mtx$v %% 1 == 0)) {
    sums  <- col_sums(mtx)
    depth <- (sum(sums) * .1) / length(sums)
    depth <- min(sums[sums >= depth])
  } else {
    depth <- 10000
  }
  
  return (as.integer(depth))
}


