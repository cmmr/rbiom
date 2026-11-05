

#' @rdname taxa_matrix
#' @export
taxa_table <- function (
    biom, rank = -1, taxa = 6, lineage = FALSE,
    md = ".all", unc = "singly", other = FALSE,
    transform = "none", ties = "random", seed = 0 ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('taxa_table', params)


  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_rank(max = Inf)
  validate_biom_field('md', max = Inf, null_ok = TRUE)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_table', params)
  tbl        <- get_cache_value(cache_file)

  if (is.null(tbl)) {

    #________________________________________________________
    # Return multiple ranks in a single table.
    #________________________________________________________
    tbl <- NULL

    for (r in rank)
      tbl %<>% dplyr::bind_rows(local({

        mtx <- taxa_matrix(
          biom      = biom,
          rank      = r,
          taxa      = taxa,
          lineage   = lineage,
          sparse    = FALSE,
          unc       = unc,
          other     = other,
          transform = transform,
          ties      = ties,
          seed      = seed )


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


    #________________________________________________________
    # Descriptive label for y-axis.
    #________________________________________________________
    resp_label <- if (eq(transform, 'percent') || biom$relative) { "Relative Abundance"
    } else if (is.null(biom$depth))                              { "Unrarefied Counts"
    } else                                                       { "Rarefied Counts"    }

    if (eq(transform, 'rank'))
      resp_label %<>% tolower() %>% paste0("Ranked Abundance\n", "(from ", ., ")")

    if (!transform %in% c('none', 'rank', 'percent'))
      resp_label %<>% paste0("\n(", transform, " transformed)")


    tbl %<>% as_rbiom_tbl()
    attr(tbl, 'response')   <- ".abundance"
    attr(tbl, 'resp_label') <- resp_label

    set_cache_value(cache_file, tbl)
  }


  return (prov_attach(tbl, cmd, biom))
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
#' @param transform   Transformation to apply to the abundance values.
#'        Options are: `c("none", "rank", "log", "log1p", "sqrt", "percent")`.
#'        `"percent"` converts counts to relative abundances before summing
#'        by taxa. `"log"` adds a pseudocount of half the smallest non-zero
#'        value in the matrix before taking the logarithm, so that absent
#'        taxa rank below the rarest observed ones. `"rank"` is useful for
#'        correcting for non-normally distributions before applying
#'        regression statistics. Default: `"none"`
#'
#' @return
#' \describe{
#'   \item{`taxa_matrix()` - }{
#'     A numeric matrix with taxa as rows, and samples as columns. }
#'   \item{`taxa_table()` - }{
#'     A tibble data frame with column names .sample, .taxa, .abundance, and any requested by `md`. }
#' }
#' Each result carries `$code`: the R code that regenerates it, including
#' the steps that produced `biom`.
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
  cmd    <- current_cmd('taxa_matrix', params)


  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_taxa(null_ok = TRUE)
  validate_rank()
  validate_bool("sparse")
  validate_var_choices('transform', c("none", "rank", "log", "log1p", "sqrt", "percent"))
  validate_var_choices('ties', c("average", "first", "last", "random", "max", "min"))
  validate_seed()


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_matrix', params)
  mtx        <- get_cache_value(cache_file)

  if (is.null(mtx)) {

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
      error = function (e) cli_abort(c("Unable to group by taxonomic level {.val {rank}}.", 'x' = conditionMessage(e))),
      expr  = local({

        counts <- if (eq(transform, 'percent')) relativize_counts(biom$counts) else biom$counts
        counts <- counts[names(map), , drop = FALSE]

        # A plain factor (no provenance class) for the S4 coercion.
        map <- if (is.factor(map)) { factor(as.character(map), levels = levels(map))
        } else                     { factor(as.character(map))                       }
        projection <- as(map, "sparseMatrix")

        mtx <- projection %*% counts
        mtx <- mtx[order(tolower(rownames(mtx))), colnames(counts), drop=FALSE]

        stopifnot(inherits(mtx, "sparseMatrix"))

        return (mtx)
      }))


    #________________________________________________________
    # Only retain taxa of interest (by abundance or name).
    #________________________________________________________
    if (!is_null(taxa)) {

      requested <- taxa
      taxa      <- top_taxa(mtx, taxa)

      if (length(taxa) == 0)
        cli_abort(c(
          'x' = "No taxa match `taxa` = {.val {requested}} at rank {.val {rank}}.",
          'i' = "Available taxa: {.val {head(rownames(mtx), 10)}}{if (nrow(mtx) > 10) ', ...' else ''}" ))


      if (isFALSE(other) || is.null(other)) {

        mtx <- mtx[taxa,,drop=FALSE]

      } else {

        if (!is_scalar_character(other) || is_na(other))
          other <- "Other"

        mtx <- mtx[setdiff(rownames(mtx), taxa),,drop=FALSE] %>%
          colSums() %>%
          matrix(., nrow = 1, dimnames = list(other, names(.))) %>%
          as("sparseMatrix") %>%
          rbind(mtx[taxa,,drop=FALSE], .)

        taxa <- c(taxa, other)
      }

    }


    #________________________________________________________
    # Optionally transform the computed abundance values.
    #________________________________________________________
    if (eq(transform, "rank")) {

      mtx %<>% as.matrix() # need those zeroes

      if (eq(ties, 'random')) {
        mtx[] <- with_seed(seed, base::rank(mtx, ties.method = ties))
      } else {
        mtx[] <- base::rank(mtx, ties.method = ties)
      }

      mtx %<>% as("sparseMatrix")

    } else if (eq(transform, "log")) {

      # Pseudocount: half the smallest non-zero value, so that zeros
      # (absent taxa) rank below the rarest observed taxa.
      mtx %<>% as.matrix()
      pseudocount <- if (any(mtx > 0)) min(mtx[mtx > 0]) / 2 else 1
      mtx <- log(mtx + pseudocount)
      mtx %<>% as("sparseMatrix")

    } else if (transform %in% c("log1p", "sqrt")) {
      mtx@x <- do.call(`::`, list('base', transform))(mtx@x)
    }


    #________________________________________________________
    # Convert to a base matrix if sparse=FALSE
    #________________________________________________________
    if (isFALSE(sparse))
      mtx <- as.matrix(mtx)

    set_cache_value(cache_file, mtx)
  }


  return (prov_result(mtx, cmd, biom))
}



#________________________________________________________
# Select taxa from a taxa x samples matrix. `taxa` is a number of taxa
# (>= 1; the most abundant), a minimum mean relative abundance (< 1), or
# a character vector of taxa names. Abundance is ranked by the mean
# relative abundance across samples. Samples with no abundance at this
# rank (e.g. from unc = "drop") are left out of the mean, rather than
# contributing 0/0 = NaN.
#________________________________________________________
top_taxa <- function (mtx, taxa) {

  if (is.character(taxa))
    return (intersect(as.character(taxa), rownames(mtx)))

  if (!is.numeric(taxa) || length(taxa) != 1)
    cli_abort("`taxa` must be a single number or a character vector, not {.type {taxa}}.")

  totals <- colSums(mtx)
  keep   <- totals > 0
  rel    <- rowMeans(t(t(mtx[, keep, drop = FALSE]) / totals[keep]))
  rel    <- sort(rel, decreasing = TRUE)

  if (taxa >= 1) { head(names(rel), taxa)
  } else         { names(rel)[rel >= taxa] }
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
  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  res    <- taxa_apply(biom, FUN = base::sum, rank, sort, lineage, unc, transform)
  prov_result(res, current_cmd('taxa_sums', params), biom)
}



#' @rdname taxa_sums
#' @export

taxa_means <- function (biom, rank = -1, sort = NULL, lineage = FALSE, unc = "singly", transform = "none") {
  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  res    <- taxa_apply(biom, FUN = base::mean, rank, sort, lineage, unc, transform)
  prov_result(res, current_cmd('taxa_means', params), biom)
}



#' @rdname taxa_sums
#' @export

taxa_apply <- function (biom, FUN, rank = -1, sort = NULL, lineage = FALSE, unc = "singly", transform = "none", ...) {

  stopifnot(is.function(FUN))
  validate_var_choices('sort', choices = c('asc', 'desc'), null_ok = TRUE)

  mtx <- taxa_matrix(biom = biom, rank = rank, lineage = lineage, sparse = TRUE, unc = unc, transform = transform)

  res <- if (identical(FUN, base::sum))  { rowSums(mtx)
  } else if (identical(FUN, base::mean)) { rowMeans(mtx)
  } else { apply(mtx, 1, FUN, ...) }

  if (is.numeric(res) && !is.null(sort))
    res <- base::sort(res, decreasing = (sort == 'desc'))

  return (res)
}





#' @rdname subset
#' @export

subset_taxa <- function (x, subset, clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) x$clone() else x
  expr <- substitute(subset)
  env  <- parent.frame()

  keep <- eval(expr = expr, envir = biom$taxonomy, enclos = env)

  if (!is.logical(keep) || length(keep) != biom$n_otus)
    cli_abort(c(
      'x' = "`subset` must evaluate to a logical vector with one value per OTU.",
      'i' = "Got {.type {keep}} of length {length(keep)}; expected length {biom$n_otus}." ))

  keep[is.na(keep)] <- FALSE

  if (!any(keep))
    cli_abort("No OTUs match `subset`.")

  pi <- prov_inline(expr, env, exclude = biom$ranks)

  biom_transform(
    biom         = biom,
    expr         = biom$counts <- biom$counts[keep, , drop = FALSE],
    call         = call('subset_taxa', quote(biom), pi$expr),
    reproducible = pi$reproducible,
    clone        = clone )
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

  keep_taxa <- split_apply(stats, '.rank', as = "list", function (x) {
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

    df <- split_apply(df, '.rank', function (x) {
      rank <- as.character(x[['.rank']][[1]])
      x[x[['.taxa']] %in% keep_taxa[[rank]], setdiff(names(x), '.rank'), drop=FALSE]
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
