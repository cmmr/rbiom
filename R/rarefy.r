
#' Rarefy OTU counts.
#' 
#' Sub-sample OTU observations such that all samples have an equal number.
#'
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' 
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
#' @param upsample   If the count data is in percentages, provide an integer 
#'        value here to scale each sample's observations to integers that sum 
#'        to this value. Generally not recommended, but can be used to 
#'        'shoehorn' metagenomic abundance estimates into rbiom's functions 
#'        that were designed for amplicon datasets. When invoked, `depth`, `n`, 
#'        and `seed` are ignored. The default, `NULL`, will throw an error if 
#'        the counts are not all integers.
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
rarefy <- function (biom, depth = 0.1, n = NULL, seed = 0L, upsample = NULL, clone = TRUE, cpus = NULL) {
  
  biom <- as_rbiom(biom)
  if (isTRUE(clone)) biom <- biom$clone()
  
  biom$counts <- mtx_rarefy(
      mtx      = biom$counts, 
      margin   = 2L,
      depth    = depth, 
      n        = n, 
      seed     = seed, 
      upsample = upsample,
      cpus     = cpus )
  
  if (isTRUE(clone)) { return (biom)            } 
  else               { return (invisible(biom)) }
}


#' Suggest a 'good' rarefaction depth.
#' 
#' @noRd
#' @keywords internal
#'     
#' @return An integer.
#'

rare_suggest <- function (mtx) {
  
  stopifnot(inherits(mtx, "sparseMatrix"))
  
  if (all(mtx@x %% 1 == 0)) {
    sums  <- colSums(mtx)
    depth <- (sum(sums) * .1) / length(sums)
    depth <- min(sums[sums >= depth])
  } else {
    depth <- 10000
  }
  
  return (as.integer(depth))
}

#' Transform a counts matrix.
#' 
#' A collection of transformations that operate directly on matrices.\cr\cr
#' Note: `rarefy_cols()`, `rescale_rows()`, and `rescale_cols()` are deprecated.
#'
#' @inherit documentation_default
#' @inherit rarefy
#' 
#' @name matrix_ops
#' @family transformations
#' @concept low_level
#'     
#' @param range   When rescaling, what should the minimum and maximum values 
#'        be? Default: `c(0, 1)`
#'
#' @param margin   Apply the transformation to the matrix's rows (`margin=1L`) 
#'        or columns (`margin=2L`). Instead of `1L` and `2L`, you may also use 
#'        `'rows'` and `'cols'`.
#'        Default: `2L` (column-wise, aka sample-wise for otu tables)
#'     
#' @return The transformed matrix. If `mtx` was a sparse matrix from the `Matrix`
#'         package, then the result will also be a sparse matrix, 
#'         otherwise the result will be a base R matrix.
#' 
#' @examples
#'     library(rbiom)
#'     
#'     # mtx_rarefy --------------------------------------
#'     biom <- hmp50$clone()
#'     sample_sums(biom) %>% head(10)
#'
#'     biom$counts %<>% mtx_rarefy(depth=1000)
#'     sample_sums(biom) %>% head(10)
#'     
#'     
#'     # rescaling ----------------------------------------
#'     mtx <- matrix(sample(1:20), nrow=4)
#'     mtx
#'     
#'     colSums(mtx)
#'     
#'     colSums(mtx_rarefy(mtx))
#'     
#'     colSums(mtx_percent(mtx))
#'     
#'     apply(mtx_rescale(mtx), 2L, max)
#'
NULL


#' @rdname matrix_ops
#' @export
mtx_rarefy <- function (mtx, margin = 2L, depth = 0.1, n = NULL, seed = 0L, upsample = NULL, cpus = NULL) {
  
  params <- eval_envir(environment())
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('mtx_rarefy', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  validate_margin()
  validate_seed()
  validate_var_range('upsample', n = 1, int = TRUE, range = c(1, Inf), null_ok = TRUE)
  validate_cpus()
  
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
  
  
  # Normally we can't rarefy non-integer counts
  if (any(mtx %% 1 != 0)) {
    
    if (is.null(upsample)) {
      
      cli_abort('Unable to rarefy non-integer counts. Set `upsample` to override.')
      
    } else {
      
      # Change percentage abundances to integer
      # Does this, but then nudges up and down until sum() == upsample:
      # round(t((t(counts) / rowSums(counts))) * upsample)
      mtx <- apply(mtx, margin, function (x) {
        
        x        <- x / sum(x)
        y        <- round(x * upsample)
        maxIters <- length(x)
        
        while (sum(y) < upsample && maxIters > 0) {
          i        <- which.min( ((y + 1) / upsample) - x )
          y[i]     <- y[i] + 1
          maxIters <- maxIters - 1
        }
        while (sum(y) > upsample && maxIters > 0) {
          i        <- which.min( x - ((y - 1) / upsample) )
          y[i]     <- y[i] - 1
          maxIters <- maxIters - 1
        }
        
        return (y)
      })
      
      if (margin == 2L) mtx <- t(mtx)
      
    }
    
  }
  
  
  # We have integer counts
  else {
    
    mtx <- ecodive::rarefy(
      counts    = mtx, 
      depth     = depth, 
      n_samples = n, 
      seed      = seed, 
      margin    = margin,
      cpus      = cpus )
  }
  
  
  set_cache_value(cache_file, mtx)
  return (mtx)
}




#' @rdname matrix_ops
#' @export
mtx_percent <- function (mtx, margin = 2L) {
  
  validate_margin()
  
  sparse <- inherits(mtx, "sparseMatrix")
  mtx    <- as.matrix(mtx)
  
  if (margin == 1) { mtx <- mtx / rowSums(mtx) }
  else             { mtx <- t(t(mtx) / colSums(mtx)) }
  
  if (sparse) mtx <- as(mtx, "sparseMatrix")
  
  return (mtx)
}




#' @rdname matrix_ops
#' @export
mtx_rescale <- function (mtx, margin = 2L, range = c(0, 1)) {
  
  validate_margin()
  validate_var_range('range',  n = 2)
  
  sparse <- inherits(mtx, "sparseMatrix")
  mtx    <- as.matrix(mtx)
  
  if (margin == 1) { mtx <- mtx / apply(mtx, 1L, max)       }
  else             { mtx <- t(t(mtx) / apply(mtx, 2L, max)) }
  
  if (!identical(range, c(0, 1))) {
    lo  <- min(range)
    hi  <- max(range)
    mtx <- mtx * (hi - lo) + lo
  }
  
  if (sparse) mtx <- as(mtx, "sparseMatrix")
  
  return (mtx)
}




#' @rdname matrix_ops
#' @export
rarefy_cols <- function (mtx, depth = 0.1, n = NULL, seed = 0L, cpus = NULL) {
  lifecycle::deprecate_soft("2.3.0", "rarefy_cols()", "mtx_rarefy()")
  mtx_rarefy(mtx = mtx, margin = 2L, depth = depth, n = n, seed = seed, cpus = cpus)
}

#' @rdname matrix_ops
#' @export
rescale_rows <- function (mtx) {
  lifecycle::deprecate_soft("2.3.0", "rescale_rows()", "mtx_percent(margin = 1L)")
  mtx_percent(mtx, margin = 1L)
}

#' @rdname matrix_ops
#' @export
rescale_cols <- function (mtx) {
  lifecycle::deprecate_soft("2.3.0", "rescale_cols()", "mtx_percent(margin = 2L)")
  mtx_percent(mtx, margin = 2L)
}
