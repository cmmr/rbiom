#' Rarefy Counts to a Constant Depth
#'
#' Normalizes the library sizes of a dataset by randomly sub-sampling
#' observations from each sample to a specific depth.
#'
#' @inherit documentation_default
#' @inherit documentation_return.biom return
#'
#' @family transformations
#' @seealso [`suggest_rarefy_depth()`] for details on the default depth selection.
#'
#' @param depth   The number of observations to keep per sample. Must be an
#'        integer greater than 0.
#'        \itemize{
#'          \item If `NULL` (the default), a depth is automatically selected
#'                that retains at least 10% of the dataset's total abundance
#'                while maximizing the number of samples kept. See
#'                [`suggest_rarefy_depth()`] for the specific heuristic used.
#'                The selected depth is recorded in `biom$code`.
#'          \item Samples with total counts less than `depth` will be
#'                dropped from the result.
#'        }
#'
#' @param seed   Random seed for selecting which observations to keep. The
#'        same `seed`, `depth`, and input always give the same result.
#'        Must be a non-negative integer. Default: `0`
#'
#' @param inflate   Logical. Handling for non-integer data (e.g. relative abundances).
#'        \itemize{
#'          \item `FALSE` (Default): The function will error if non-integers are detected.
#'                Rarefaction requires discrete counts (integers).
#'          \item `TRUE`: The function will automatically rescale (inflate)
#'                non-integers to integers using [`biom_inflate()`] before rarefying.
#'                This is useful for 'shoehorning' metagenomic relative abundance
#'                data into diversity functions that strictly require integers.
#'        }
#'
#' @description
#' This function reduces the number of observations (reads) in each sample to
#' a fixed integer value (`depth`). Samples with fewer observations than the
#' specified depth are discarded.
#'
#' Rarefaction is a common technique in microbiome analysis used to account for
#' uneven sequencing effort across samples. By standardizing the library size,
#' it allows for fair comparisons of alpha and beta diversity metrics.
#'
#' The transformation is recorded in the object's history (`biom$code`).
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- hmp50[1:5]
#'     sample_sums(biom)
#'
#'     # Rarefy to the lowest sample depth
#'     # (All samples are kept, but counts are reduced)
#'     biom_rare <- rarefy(biom, depth = min(sample_sums(biom)))
#'     sample_sums(biom_rare)
#'
#'     # Auto-select depth (may drop samples with low coverage)
#'     biom_auto <- rarefy(biom)
#'     sample_sums(biom_auto)
#'
#'     biom_auto$code
#'
rarefy <- function(
    biom, depth = NULL, seed = 0L, inflate = FALSE, clone = TRUE, cpus = n_cpus()) {

  biom <- as_rbiom(biom)
  if (isTRUE(clone)) biom <- biom$clone()

  validate_var_range('depth', int = TRUE, range = c(1, Inf), null_ok = TRUE, n = 1)
  validate_seed()
  validate_bool('inflate')
  validate_cpus()


  #________________________________________________________
  # Non-integer counts must be inflated first.
  #________________________________________________________
  counts  <- biom$counts
  non_int <- any(counts@x %% 1 != 0)

  if (non_int) {
    if (!isTRUE(inflate))
      cli_abort(c(
        'x' = "Unable to rarefy non-integer counts.",
        'i' = "Set `inflate = TRUE` to scale them to integers first (see `biom_inflate()`)." ))
    counts <- inflate_counts(counts)
  }


  #________________________________________________________
  # Default depth selection.
  #________________________________________________________
  if (is.null(depth))
    depth <- depth_from_sums(colSums(counts))

  depth <- as.integer(depth)


  #________________________________________________________
  # Execute rarefaction and record it in the object's history.
  #________________________________________________________
  pc <- prov_call('rarefy', c(
    list(quote(biom), depth = depth),
    if (!eq(seed, 0L))     list(seed    = seed),
    if (isTRUE(inflate))   list(inflate = TRUE) ))

  biom_transform(
    biom = biom,
    expr = biom$counts <- ecodive::rarefy(
      counts = counts,
      margin = 2L,
      depth  = depth,
      seed   = seed,
      cpus   = cpus,
      drop   = TRUE,
      warn   = FALSE ),
    call         = pc$call,
    reproducible = pc$reproducible,
    clone        = clone )
}


#' Inflate Relative Abundances to Integer Counts
#'
#' Scaling a matrix of proportions (or counts) to a new target depth,
#' rounding to integers while preserving the original total abundance sum
#' exactly.
#'
#' @inherit documentation_default
#' @inherit documentation_return.biom return
#'
#' @family transformations
#' @seealso [`suggest_inflate_depths()`] for details on how target depths are
#'          estimated when `depth = NULL`.
#'
#' @param depth   The target library size (sum) for each sample: a single
#'        integer, or one integer per sample. Must be greater than 0. If
#'        `NULL` (the default), the depth is estimated per-sample using the
#'        "Singleton Peak Heuristic". See [`suggest_inflate_depths()`] for
#'        algorithm details. Counts that are already integers are left
#'        unchanged when `depth = NULL`.
#'
#' @section Rounding (Largest Remainder Method):
#'
#' To ensure the sum of the resulting counts equals the target `depth` exactly
#' (avoiding drift caused by simple rounding), this function uses the Largest
#' Remainder Method (also known as the Hare-Niemeyer method).
#'
#' It assigns the integer part of the scaled value to each feature, and then
#' distributes the remaining counts to the features with the largest fractional
#' parts.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- hmp50[1:5]
#'     sample_sums(biom)
#'
#'     biom <- biom_relativize(biom)
#'     sample_sums(biom)
#'
#'     biom <- biom_inflate(biom)
#'     sample_sums(biom)
#'
biom_inflate <- function(biom, depth = NULL, clone = TRUE) {

  biom <- as_rbiom(biom)
  if (isTRUE(clone)) biom <- biom$clone()

  validate_var_range('depth', int = TRUE, range = c(1, Inf), null_ok = TRUE)

  if (!is.null(depth) && !length(depth) %in% c(1L, biom$n_samples))
    cli_abort("`depth` must be a single value or one value per sample ({biom$n_samples}), not {length(depth)}.")

  pc <- prov_call('biom_inflate', c(list(quote(biom)), if (!is.null(depth)) list(depth = depth)))

  biom_transform(
    biom         = biom,
    expr         = biom$counts <- inflate_counts(biom$counts, depth),
    call         = pc$call,
    reproducible = pc$reproducible,
    clone        = clone )
}



#________________________________________________________
# Largest Remainder Method, column-wise on a dgCMatrix. `depth` is NULL
# (estimate per sample; integer counts are returned unchanged), a single
# value, or one value per column.
#________________________________________________________
inflate_counts <- function (mtx, depth = NULL) {

  if (is.null(depth)) {
    if (!any(mtx@x %% 1 != 0)) return (mtx)
    depth <- inflate_depths(mtx)
  }

  if (length(depth) == 1) depth <- rep(depth, ncol(mtx))
  stopifnot(length(depth) == ncol(mtx))

  for (i in seq_len(ncol(mtx))) {

    start <- mtx@p[i] + 1
    end   <- mtx@p[i+1]
    if (end < start) next # Empty column

    x      <- mtx@x[start:end]
    p      <- x / sum(x)
    quotas <- p * depth[[i]]
    counts <- floor(quotas)

    remainder <- depth[[i]] - sum(counts)
    if (remainder > 0) {
      fractions   <- quotas - counts
      add_indices <- tail(order(fractions), remainder)
      counts[add_indices] <- counts[add_indices] + 1
    }

    mtx@x[start:end] <- counts
  }

  return (Matrix::drop0(mtx))
}



#' Relativize Counts to Proportions
#'
#' Convert absolute counts to relative abundances (proportions) where each
#' sample sums to 1.
#'
#' @inherit documentation_default
#' @inherit documentation_return.biom return
#'
#' @family transformations
#'
#' @description
#' This function normalizes the data by dividing each observation by the total
#' library size of its sample. The resulting values represent the proportion
#' (0 to 1) of the sample composed of that specific feature.
#'
#' This is a common transformation for microbiome data, as it accounts for
#' differences in sequencing depth across samples, allowing for comparison of
#' community composition. The result has `biom$relative == TRUE`.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- hmp50[1:5]
#'
#'     # Raw counts sum to different library sizes
#'     sample_sums(biom)
#'
#'     # Relativized counts sum to 1
#'     biom_rel <- biom_relativize(biom)
#'     sample_sums(biom_rel)
#'     biom_rel$relative
#'
biom_relativize <- function(biom, clone = TRUE) {

  biom <- as_rbiom(biom)
  if (isTRUE(clone)) biom <- biom$clone()

  biom_transform(
    biom  = biom,
    expr  = biom$counts <- relativize_counts(biom$counts),
    call  = quote(biom_relativize(biom)),
    clone = clone )
}



#________________________________________________________
# Column-normalize a dgCMatrix so that each column sums to 1.
#________________________________________________________
relativize_counts <- function (mtx) {

  sums <- Matrix::colSums(mtx)
  sums[sums == 0] <- 1

  res <- mtx %*% Matrix::Diagonal(x = 1/sums)
  dimnames(res) <- dimnames(mtx)

  return (as(res, 'dgCMatrix'))
}



#' Rescale Counts to a Specific Range
#'
#' Linearly rescale each sample's values to lie between a specified minimum
#' and maximum.
#'
#' @inherit documentation_default
#' @inherit documentation_return.biom return
#'
#' @family transformations
#'
#' @param range   Numeric vector of length 2. Target min and max.
#'        Default: `c(0, 1)`.
#'
#' @description
#' This function performs a min-max scaling on each sample independently.
#'
#' It is useful for normalization techniques that require data to be within a
#' specific bounded range, or for visualization purposes where maintaining the
#' relative distances between values is important but the absolute magnitude
#' needs adjustment.
#'
#' @section Mathematical Transformation:
#'
#' The rescaling is performed in two steps:
#'
#' 1.  **Normalize:** Divide values by the maximum value in that sample,
#'     scaling them to a `[0, 1]` range relative to the sample's peak.
#'
#' 2.  **Scale and Shift:** Apply the target range using the formula:
#'     \deqn{x_{new} = x_{norm} \times (max - min) + min}
#'
#' @note
#' If `range` starts at a non-zero value (e.g., `c(1, 10)`), the sparsity of the
#' matrix will be destroyed because all zero counts will be shifted to the
#' minimum value. This can significantly increase memory usage for large datasets.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- hmp50[1:5]
#'
#'     # Original range
#'     range(as.matrix(biom))
#'
#'     # Rescaled to 0-1
#'     biom_01 <- biom_rescale(biom)
#'     range(as.matrix(biom_01))
#'
#'     # Rescaled to 0-100 (Percentages)
#'     biom_100 <- biom_rescale(biom, range = c(0, 100))
#'     range(as.matrix(biom_100))
#'
biom_rescale <- function(biom, range = c(0, 1), clone = TRUE) {

  biom <- as_rbiom(biom)
  if (isTRUE(clone)) biom <- biom$clone()

  validate_var_range('range', n = 2)

  pc <- prov_call('biom_rescale', c(list(quote(biom)), if (!identical(range, c(0, 1))) list(range = range)))

  biom_transform(
    biom = biom,
    expr = {

      # 1. Scale to 0-1 (Divide by Max per sample/column)
      maxs <- apply(biom$counts, 2, max)
      maxs[maxs == 0] <- 1
      scl_counts <- biom$counts %*% Matrix::Diagonal(x = 1/maxs)

      # 2. Rescale to target range
      if (!identical(range, c(0, 1))) {
        lo  <- min(range)
        hi  <- max(range)
        # Note: Adding 'lo' destroys sparsity if lo != 0
        scl_counts <- scl_counts * (hi - lo) + lo
      }

      dimnames(scl_counts) <- dimnames(biom$counts)
      biom$counts <- scl_counts
    },
    call         = pc$call,
    reproducible = pc$reproducible,
    clone        = clone )
}


#' Suggest Rarefaction Depth
#'
#' Calculates a rarefaction depth that balances retaining samples against
#' retaining total observations.
#'
#' @inherit documentation_default
#'
#' @seealso [`rarefy()`] which uses this heuristic when `depth = NULL`.
#'
#' @section Heuristic:
#' This function selects a depth by analyzing the trade-off between dropping
#' samples (to increase depth) and lowering depth (to keep samples).
#'
#' 1.  **Calculate Yields:** For every distinct sample depth in the dataset,
#'     calculate the total number of observations that would remain if the
#'     dataset were rarefied to that level.
#'     \deqn{Yield_d = d \times N_{\ge d}}
#'     Where \eqn{d} is the depth and \eqn{N_{\ge d}} is the number of samples
#'     with at least that many reads.
#'
#' 2.  **Define Threshold:** Calculate 10% of the total observations in the
#'     original un-rarefied dataset.
#'
#' 3.  **Select Depth:** Find the **lowest** depth \eqn{d} where the \eqn{Yield_d}
#'     exceeds this 10% threshold.
#'
#' This approach prioritizes keeping as many samples as possible, provided that
#' doing so doesn't discard more than 90% of the dataset's total information.
#'
#' @return A single integer representing the suggested rarefaction depth.
#' @export
#' @examples
#'     library(rbiom)
#'
#'     suggest_rarefy_depth(hmp50)
#'
suggest_rarefy_depth <- function (biom) {

  biom <- as_rbiom(biom)
  depth_from_sums(Matrix::colSums(biom$counts))
}



#________________________________________________________
# The rarefaction heuristic, from a vector of sample sums.
#________________________________________________________
depth_from_sums <- function (sums) {

  sums_sorted <- sort(sums)
  total_reads <- sum(sums_sorted)
  if (total_reads == 0)
    cli_abort("Matrix contains no observations.")

  # Select the lowest depth that still retains at least 10%
  # of the total observations in the dataset.
  n_sams    <- length(sums_sorted)
  threshold <- 0.1 * total_reads

  # Find the first index (lowest depth) that meets the threshold
  yields <- sums_sorted * seq.int(n_sams, 1L)
  idx    <- which(yields >= threshold)[1L]

  if (!is.na(idx)) { depth <- sums_sorted[idx] }
  else             { depth <- max(sums_sorted) } # nocov

  return (as.integer(depth))
}


#' Suggest Inflation Depths
#'
#' Estimates the optimal sequencing depth for each sample in a matrix by
#' leveraging the global abundance distribution structure.
#'
#' @inherit documentation_default
#'
#' @seealso [`biom_inflate()`] which uses this heuristic when `depth = NULL`.
#'
#' @param adjust   Numeric. Bandwidth adjustment for the kernel density
#'        estimation. Default: `1.5`.
#'
#' @section The Singleton Peak Heuristic:
#'
#' When `depth = NULL`, [`biom_inflate()`] calls this function to estimate the original
#' sequencing depth for each sample. The underlying assumption is that in
#' typical microbiome datasets, the most frequent count value (the mode of the
#' abundance distribution) is 1 (a singleton).
#'
#' The algorithm works as follows:
#'
#' 1.  **Log-Transformation:** Non-zero relative abundances are log10-transformed.
#'
#' 2.  **Global Consensus:** To overcome sparsity in individual samples, distributions
#'     are centered by their medians and aggregated across all samples.
#'
#' 3.  **Peak Detection:** Kernel Density Estimation (KDE) is used to identify the
#'     peak (mode) of this aggregated distribution.
#'
#' 4.  **Scaling:** A scaling factor is calculated for each sample that shifts
#'     this peak to correspond to an integer count of 1.
#'
#' This approach effectively "shoehorns" relative abundance data into integer
#' formats required by diversity metrics (like rarefaction or Chao1) by
#' maximizing the number of singletons in the resulting matrix.
#'
#' Samples with fewer than two non-zero features carry no information about
#' the distribution's shape; they are assigned the median depth of the other
#' samples (or, if there are too few informative samples for the heuristic,
#' every sample is scaled so that its smallest non-zero value becomes 1).
#' Counts that are already integers are returned as their current sums.
#'
#' @return A named integer vector of recommended depths for each sample.
#' @export
#' @examples
#'     library(rbiom)
#'
#'     depths <- suggest_inflate_depths(biom_relativize(hmp50))
#'     head(depths)
#'
suggest_inflate_depths <- function(biom, adjust = 1.5) {

  biom <- as_rbiom(biom)
  mtx  <- biom$counts

  if (!any(mtx@x %% 1 != 0))
    return (stats::setNames(as.integer(round(Matrix::colSums(mtx))), colnames(mtx)))

  inflate_depths(mtx, adjust)
}



#________________________________________________________
# The singleton peak heuristic on a dgCMatrix.
#________________________________________________________
inflate_depths <- function (mtx, adjust = 1.5) {

  n_sams <- ncol(mtx)
  sams   <- colnames(mtx)

  global_centered_vals <- numeric(0)
  sample_medians       <- numeric(n_sams)
  sample_minimums      <- rep(NA_real_, n_sams)
  valid_samples        <- logical(n_sams)


  # Iterate samples using the sparse structure: column i's values are
  # x[(p[i] + 1):p[i + 1]].
  p <- mtx@p
  v <- mtx@x

  for (i in seq_len(n_sams)) {

    if (p[i+1] == p[i]) next # Empty column

    vals <- v[ (p[i] + 1) : p[i+1] ]
    vals <- vals[is.finite(vals) & vals > 0]
    if (length(vals) == 0) next

    sample_minimums[i] <- min(vals)
    if (length(vals) < 2) next

    lvals <- log10(vals)
    med   <- median(lvals)
    sample_medians[i]    <- med
    valid_samples[i]     <- TRUE
    global_centered_vals <- c(global_centered_vals, lvals - med)
  }


  results <- rep(NA_integer_, n_sams)
  if (!is.null(sams)) names(results) <- sams


  if (length(global_centered_vals) >= 10) {

    # Global Consensus KDE
    bw <- tryCatch(stats::bw.SJ(global_centered_vals), error = function(e) "nrd0")
    bw_val <- if(is.numeric(bw)) bw else stats::bw.nrd0(global_centered_vals)

    step_size <- (bw_val * adjust) / 10
    span      <- diff(range(global_centered_vals))
    n_grid    <- min(4096, max(512, 2^ceiling(log2(span / step_size))))

    dens <- stats::density(global_centered_vals, bw = bw, adjust = adjust, n = n_grid)

    # Find Peak
    y_diff       <- diff(sign(diff(dens$y)))
    peak_indices <- which(y_diff == -2) + 1
    if (length(peak_indices) == 0) peak_indices <- which.max(dens$y)

    peak_indices     <- peak_indices[dens$y[peak_indices] > (0.10 * max(dens$y))]
    best_idx         <- peak_indices[which.min(dens$x[peak_indices])]
    singleton_offset <- dens$x[best_idx]

    # Project
    for (i in which(valid_samples)) {
      est_peak   <- sample_medians[i] + singleton_offset
      results[i] <- as.integer(round(1 / 10^est_peak))
    }
  }


  # Samples the heuristic can't handle: use the median of the others,
  # or make the smallest observed value a singleton.
  fallback <- is.na(results) & !is.na(sample_minimums)

  if (any(fallback)) {
    if (any(!is.na(results))) {
      results[fallback] <- as.integer(round(median(results, na.rm = TRUE)))
    } else {
      results[fallback] <- as.integer(round(1 / sample_minimums[fallback]))
    }
  }

  results[!is.na(results) & results < 1L] <- 1L

  return (results)
}





# Deprecated ------------------------------------------------------------------

#' Deprecated matrix transformations
#'
#' A collection of transformations that operate directly on matrices.
#'
#' @inherit documentation_default
#'
#' @name matrix_ops
#'
#' @param mtx   A numeric matrix or sparse matrix of counts.
#'
#' @param margin   Apply the transformation to the matrix's rows (`margin=1L`)
#'        or columns (`margin=2L`). Instead of `1L` and `2L`, you may also use
#'        `'rows'` and `'cols'`.
#'        Default: `2L` (column-wise, aka sample-wise for otu tables)
#'
#' @param range   When rescaling, what should the minimum and maximum values
#'        be? Default: `c(0, 1)`
#'
#' @param depth   How many observations to keep per sample.
#'
#' @param n   Deprecated. The number of samples to keep. This argument is
#'        ignored in the current version.
#'
#' @param seed   An integer seed for randomizing which observations to keep or
#'        drop.
#'
#' @param upsample   If the count data is in percentages, provide an integer
#'        value here to scale each sample's observations to integers that sum
#'        to this value before rarefying. See [biom_inflate()].
#'
#'
#' @return A A numeric matrix or sparse matrix, depending on the input type,
#'         with the same dimensions as `mtx`.
#'
#' @examples
#'     mtx <- matrix(1:12, ncol = 3, dimnames = list(paste0("OTU", 1:4), paste0("Sample", 1:3)))
#'
#'     mtx
#'
#'     suppressWarnings({
#'
#'       mtx_rarefy(mtx)
#'
#'       rarefy_cols(mtx)
#'
#'       mtx_percent(mtx)
#'
#'       mtx_rescale(mtx)
#'
#'       rescale_rows(mtx)
#'
#'       rescale_cols(mtx)
#'
#'     })
#'
NULL


#' @rdname matrix_ops
#' @export
mtx_rarefy <- function(
    mtx, margin = 2L, depth = 0.1, n = NULL, seed = 0L, upsample = NULL, cpus = NULL) {

  lifecycle::deprecate_soft("3.0.0", "mtx_rarefy()", "rarefy()")

  # 1. Handle Orientation (New functions expect samples as columns)
  # If margin=1 (rows are samples), we must transpose inputs.
  is_row_margin <- (margin == 1L || margin == 'rows')
  if (is_row_margin) mtx <- t(mtx)

  biom <- as_rbiom(mtx)

  # 2. `upsample` scales each sample to this many integer counts first.
  if (!is.null(upsample))
    biom <- biom_inflate(biom, depth = upsample, clone = FALSE)

  # 3. Handle Depth Logic (Old allowed < 1 for auto-calc, New requires int >= 1)
  # If depth is a fraction, calculate a specific integer depth based on the input
  if (is.numeric(depth) && length(depth) == 1 && depth > 0 && depth < 1) {
    sums  <- Matrix::colSums(biom$counts)
    depth <- as.integer((sum(sums) * depth) / length(sums))
    if (depth < 1) depth <- 1L
  }

  # 4. Call New Implementation
  # Note: 'n' is effectively ignored as the new implementation relies on 'depth'
  # to determine sample retention (via ecodive::rarefy's internal logic).
  biom <- rarefy(
    biom    = biom,
    depth   = depth,
    seed    = seed,
    inflate = FALSE,
    clone   = FALSE,
    cpus    = if (is.null(cpus)) n_cpus() else cpus )

  # 5. Return to original format
  res <- biom$counts
  if (is_row_margin) res <- t(res)

  # Convert to base matrix if input wasn't sparse, to match old behavior often expecting matrix
  if (!inherits(mtx, "sparseMatrix"))
    res <- as.matrix(res)

  return(res)
}


#' @rdname matrix_ops
#' @export
mtx_percent <- function(mtx, margin = 2L) {

  lifecycle::deprecate_soft("3.0.0", "mtx_percent()", "biom_relativize()")

  # Handle Orientation
  is_row_margin <- (margin == 1L || margin == 'rows')
  if (is_row_margin) mtx <- t(mtx)

  # Execute
  biom <- biom_relativize(biom = mtx, clone = TRUE)

  # Restore Orientation
  res <- biom$counts
  if (is_row_margin) res <- t(res)

  # Format Match
  if (!inherits(mtx, "sparseMatrix"))
    res <- as.matrix(res)

  return(res)
}


#' @rdname matrix_ops
#' @export
mtx_rescale <- function(mtx, margin = 2L, range = c(0, 1)) {

  lifecycle::deprecate_soft("3.0.0", "mtx_rescale()", "biom_rescale()")

  # Handle Orientation
  is_row_margin <- (margin == 1L || margin == 'rows')
  if (is_row_margin) mtx <- t(mtx)

  # Execute
  biom <- biom_rescale(biom = mtx, range = range, clone = TRUE)

  # Restore Orientation
  res <- biom$counts
  if (is_row_margin) res <- t(res)

  # Format Match
  if (!inherits(mtx, "sparseMatrix"))
    res <- as.matrix(res)

  return(res)
}


#' @rdname matrix_ops
#' @export
rarefy_cols <- function(mtx, depth = 0.1, n = NULL, seed = 0L, cpus = NULL) {
  lifecycle::deprecate_soft("3.0.0", "rarefy_cols()", "rarefy()")
  mtx_rarefy(mtx = mtx, margin = 2L, depth = depth, n = n, seed = seed, cpus = cpus)
}

#' @rdname matrix_ops
#' @export
rescale_rows <- function(mtx) {
  lifecycle::deprecate_soft("3.0.0", "rescale_rows()", "biom_rescale()")
  mtx_percent(mtx, margin = 1L)
}

#' @rdname matrix_ops
#' @export
rescale_cols <- function(mtx) {
  lifecycle::deprecate_soft("3.0.0", "rescale_cols()", "biom_rescale()")
  mtx_percent(mtx, margin = 2L)
}
