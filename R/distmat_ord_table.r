
#' Run ordinations on a distance matrix.
#'
#' @inherit documentation_default
#'
#' @family ordination
#'
#' @param ord    Method for reducing dimensionality. Options are:
#'        \describe{
#'            \item{`"PCoA"` - }{ Principal coordinate analysis; [stats::cmdscale()]. }
#'            \item{`"UMAP"` - }{ Uniform manifold approximation and projection; [uwot::umap()]. }
#'            \item{`"NMDS"` - }{ Nonmetric multidimensional scaling; [vegan::metaMDS()]. }
#'            \item{`"tSNE"` - }{ t-distributed stochastic neighbor embedding; [Rtsne::Rtsne()]. }
#'        }
#'        Multiple/abbreviated values allowed. Default: `"PCoA"`
#'
#' @param seed   Random seed for the stochastic ordinations (NMDS, tSNE,
#'        UMAP), so that they are reproducible. PCoA is deterministic.
#'        Must be a non-negative integer. Default: `0`
#'
#' @param ...  Additional arguments for the `ord` function. For tSNE,
#'        `perplexity` defaults to 30, or to `(n - 1) / 3` for fewer than
#'        91 samples.
#'
#' @return A tibble data.frame with columns \code{.ord}, \code{.sample},
#'         \code{.x}, \code{.y}, and (optionally) \code{.z}. The attribute
#'         `'axes'` is a tibble with one row per axis (`.ord`, `.axis`,
#'         `.label`, `.var`): for PCoA, `.var` is the percent of the
#'         variance (the positive eigenvalues) explained by that axis and
#'         `.label` reads e.g. `"PCoA 1 (24.3%)"`; for the other ordinations
#'         `.var` is `NA`.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     dm  <- bdiv_distmat(hmp50, "bray")
#'     ord <- distmat_ord_table(dm, "PCoA")
#'     head(ord)
#'
#'     attr(ord, 'axes')
#'

distmat_ord_table <- function (dm, ord = "PCoA", k = 2L, seed = 0, ...) {

  params <- eval_envir(environment(), ...)

  cmd <- local({
    p <- as.list(params)
    p$dm <- aa(NA, display = "dm")
    current_cmd('distmat_ord_table', p)
  })
  schema_check('distmat_ord_table', params)


  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  if (!inherits(dm, 'dist'))
    cli_abort("`dm` must be a dist object, not {.type {dm}}.")

  if (attr(dm, 'Size') < 3)
    cli_abort("At least three samples are needed for an ordination; `dm` has {attr(dm, 'Size')}.")

  validate_ord(max = Inf)
  validate_var_range('k', range = c(2, 3), n = 1, int = TRUE)
  validate_seed()


  #________________________________________________________
  # See if this result is already in the cache. Provenance attributes
  # on `dm` are excluded from the key so that identical distances from
  # different pipelines share an entry.
  #________________________________________________________
  cache_file <- local({
    key <- as.list(params, all.names = TRUE)
    key$dm <- strip_prov(dm)
    get_cache_file('distmat_ord_table', key)
  })
  tbl        <- get_cache_value(cache_file)

  if (is.null(tbl)) {

    dots   <- params$.dots
    labels <- attr(dm, "Labels", exact = TRUE)
    axes   <- list()


    #________________________________________________________
    # Call 3rd party function with proper arguments. Stochastic
    # methods run with the RNG seeded by `seed`.
    #________________________________________________________
    tbl <- dplyr::bind_rows(lapply(ord, function (o) {

      var <- rep(NA_real_, k)

      if (o == "PCoA") {
        args <- c(list(d = dm, k = k, eig = TRUE), fun_params(stats::cmdscale, dots)[setdiff(names(dots), c('d', 'k', 'eig'))])
        pc   <- do.call(stats::cmdscale, args)
        res  <- pc$points
        eig  <- pc$eig
        var  <- 100 * eig[seq_len(k)] / sum(eig[eig > 0])

      } else if (o == "tSNE") {
        # Rtsne requires perplexity <= (n - 1) / 3.
        require_package('Rtsne', 'to use "tSNE" ordination')
        args <- c(list(X = dm, dims = k), fun_params(Rtsne::Rtsne, dots)[setdiff(names(dots), c('X', 'dims', 'is_distance'))])
        args[['perplexity']] %<>% if.null(min(30, (attr(dm, 'Size') - 1) / 3))
        res  <- with_seed(seed, do.call(Rtsne::Rtsne, args))[['Y']]
        rownames(res) <- labels

      } else if (o == "NMDS") {
        args <- c(list(comm = dm, k = k), fun_params(vegan::metaMDS, dots)[setdiff(names(dots), c('comm', 'k'))])
        args[['trace']] %<>% if.null(0)
        res  <- with_seed(seed, do.call(vegan::metaMDS, args))[['points']]
        rownames(res) <- labels

      } else if (o == "UMAP") {
        require_package('uwot', 'to use "UMAP" ordination')
        args <- c(list(X = dm), fun_params(uwot::umap, dots)[setdiff(names(dots), 'X')])
        args[['n_components']] %<>% if.null(k)
        args[['n_neighbors']]  %<>% if.null(max(2, min(100, as.integer(attr(dm, 'Size') / 3))))
        res  <- with_seed(seed, do.call(uwot::umap, args))
        rownames(res) <- labels
      }

      axis_names <- head(c(".x", ".y", ".z"), k)
      axes[[o]] <<- tibble(
        .ord   = o,
        .axis  = axis_names,
        .label = ifelse(is.na(var), paste(o, seq_len(k)), sprintf("%s %d (%.1f%%)", o, seq_len(k), var)),
        .var   = var )

      res <- as.data.frame(res[, seq_len(k), drop = FALSE])
      colnames(res) <- axis_names

      tibble(.ord = o, .sample = rownames(res), res)
    }))

    tbl[['.ord']] %<>% factor(levels = ord)

    axes <- dplyr::bind_rows(axes)
    axes[['.ord']] %<>% factor(levels = ord)


    #________________________________________________________
    # Table header.
    #________________________________________________________
    tbl %<>% as_rbiom_tbl()
    attr(tbl, 'axes') <- axes

    if (length(ord) == 1)
      attr(tbl, 'tbl_sum') <- c('Ordination' = switch(
        EXPR = ord,
        PCoA = sprintf("stats::cmdscale(dm, k = %s, eig = TRUE)", k),
        tSNE = sprintf("Rtsne::Rtsne(dm, dims = %s)", k),
        NMDS = sprintf("vegan::metaMDS(dm, k = %s)", k),
        UMAP = sprintf("uwot::umap(dm, n_components = %s)", k) ))

    set_cache_value(cache_file, tbl)
  }


  return (prov_attach(tbl, cmd, data = dm))
}
