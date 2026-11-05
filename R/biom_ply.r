#' Apply a function to each subset of an rbiom object.
#'
#' `blply()` and `bdply()` let you divide your biom dataset into smaller
#' pieces, run a function on those smaller rbiom objects, and return the
#' results as a data.frame or list.
#'
#' You can also specify additional variables for your function to iterate over
#' in unique combinations.
#'
#' Each subset passed to `FUN` is a proper rbiom object whose history records
#' how it was made, e.g. `subset(biom, Sex == "Male")`, so that results
#' computed from it (and their `$code`) are reproducible.
#'
#' @inherit documentation_default
#'
#' @family metadata
#' @family biom
#'
#' @param vars   A character vector of metadata fields. Each unique combination
#'        of values in these columns will be used to create a subsetted
#'        rbiom object to pass to `FUN.` If `NULL`,
#'        `biom` will be passed to `FUN` unaltered. Unambiguous
#'        abbreviations of metadata fields are also accepted. Samples with an
#'        `NA` in any of these fields are left out.
#'
#' @param FUN   The function to execute on each subset of `biom`.
#'        For `bdply()`, the returned value will be coerced to a data.frame.
#'        For `blply()`, any returned value is unmodified.
#'
#' @param ...   Additional arguments to pass on to `FUN`.
#'
#' @param iters   A named list of values to pass to `FUN`. Unlike
#'        \code{...}, these will be iterated over in all combinations.
#'        Default: \code{list()}
#'
#' @param prefix   When `TRUE`, prefixes the names in in \code{iters} with
#'        a '.' in the final data.frame or 'split_labels' attribute.
#'        Default: `FALSE`
#'
#' @return For `bdply()`, a tibble data.frame comprising the accumulated
#'         outputs of `FUN`, along with the columns specified by
#'         `vars` and `iters`. For `blply()`, a named list that has details
#'         about `vars` and `iters` in `attr(,'split_labels')`.
#'
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     bdply(hmp50, "Sex", `$`, 'n_samples')
#'
#'     blply(hmp50, "Sex", `$`, 'n_samples') %>% unlist()
#'
#'     bdply(hmp50, c("Body Site", "Sex"), function (b) {
#'       adm <- adiv_matrix(b)[,c("shannon", "simpson")]
#'       apply(adm, 2L, mean)
#'     })
#'
#'     iters <- list(d = c("bray", "euclid"))
#'     bdply(hmp50, "Sex", iters = iters, function (b, d) {
#'       r <- range(bdiv_distmat(biom = b, bdiv = d))
#'       round(data.frame(min = r[[1]], max = r[[2]]))
#'     })
#'
bdply <- function (biom, vars, FUN, ..., iters = list(), prefix = FALSE) {


  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  if (!inherits(FUN, 'function'))
    cli_abort("`FUN` must be a function, not {.type {FUN}}.")

  dots  <- list(...)
  iters <- expand.grid(iters, stringsAsFactors = FALSE)

  apply_iters <- function (b) {
    if (nrow(iters) == 0) return (do.call(FUN, c(list(b), dots)))
    rows_apply(iters, function (iter) do.call(FUN, c(list(b), dots, as.list(iter))))
  }


  #________________________________________________________
  # Simple cases where we're not faceting by metadata.
  #________________________________________________________
  if (is.null(vars)) {

    result <- apply_iters(biom)
    if (!is.data.frame(result)) result <- as_tibble(result)

  } else {

    biom <- ply_biom(biom)
    validate_biom_field('vars', col_type = 'cat', max = Inf, null_ok = TRUE)

    result <- split_apply(biom$metadata, vars, function (df) {
      apply_iters(ply_subset(biom, df, vars))
    })
  }


  if (isTRUE(prefix) && nrow(iters) > 0) {
    i <- seq_len(ncol(iters)) + length(vars)
    colnames(result)[i] <- paste0(".", colnames(iters))
  }


  return (as_rbiom_tbl(result))
}



#' @rdname bdply
#' @export

blply <- function (biom, vars, FUN, ..., iters = list(), prefix = FALSE) {


  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  if (!inherits(FUN, 'function'))
    cli_abort("`FUN` must be a function, not {.type {FUN}}.")

  dots  <- list(...)
  iters <- expand.grid(iters, stringsAsFactors = FALSE)

  apply_iters <- function (b) {
    if (nrow(iters) == 0) return (do.call(FUN, c(list(b), dots)))
    res <- lapply(seq_len(nrow(iters)), function (i) {
      do.call(FUN, c(list(b), dots, as.list(iters[i, , drop = FALSE]))) })
    names(res) <- seq_len(nrow(iters))
    return (res)
  }


  #________________________________________________________
  # Simple cases where we're not faceting by metadata.
  #________________________________________________________
  if (is.null(vars)) {

    if (nrow(iters) == 0) {
      result <- list(apply_iters(biom))
    } else {
      result <- apply_iters(biom)
      attr(result, 'split_labels') <- iters
    }

  } else {

    biom <- ply_biom(biom)
    validate_biom_field('vars', col_type = 'cat', max = Inf, null_ok = TRUE)

    result <- split_apply(biom$metadata, vars, as = "list", function (df) {
      apply_iters(ply_subset(biom, df, vars))
    })


    #________________________________________________________
    # Un-nest the per-stratum lists of iterations.
    #________________________________________________________
    if (nrow(iters) > 0) {

      split_labels <- attr(result, 'split_labels', exact = TRUE)
      result       <- do.call(c, result)

      attr(result, 'split_labels') <- rows_apply(split_labels, function (row) iters)
    }
  }


  if (isTRUE(prefix) && nrow(iters) > 0)
    attr(result, 'split_labels') <- local({
      df <- attr(result, 'split_labels', exact = TRUE)
      i  <- seq_len(ncol(iters)) + length(vars)
      colnames(df)[i] <- paste0(".", colnames(iters))
      return (df)
    })


  return (result)
}



#________________________________________________________
# Code text that assigns `stats` from `inner` (lines of code whose last
# top-level expression is the value), run for each `split.by` stratum and
# `bdiv` metric with bdply(), or directly when there is only one of each.
# `inner` may refer to `biom` and, when several metrics are given, `bdiv`.
#________________________________________________________
bdply_code <- function (inner, split.by, bdiv, var = "stats", iterate = length(bdiv) > 1) {

  if (!is.null(split.by) || isTRUE(iterate)) {

    args <- c(
      "biom",
      sprintf("vars = %s", if (is.null(split.by)) "NULL" else as.args(list(split.by))),
      if (isTRUE(iterate)) sprintf("iters = list(bdiv = %s), prefix = TRUE", as.args(list(bdiv))),
      sprintf("FUN = function (%s) {", if (isTRUE(iterate)) "biom, bdiv" else "biom") )

    c(sprintf("%s <- bdply(%s", var, paste(args, collapse = ", ")), paste0("  ", inner), "})")

  } else {

    i <- max(which(!startsWith(inner, " ") & !startsWith(inner, "}")))
    inner[[i]] <- sprintf("%s <- %s", var, inner[[i]])
    inner
  }
}



#________________________________________________________
# Coerce to an rbiom object, with a clear error if that's impossible.
#________________________________________________________
ply_biom <- function (biom) {

  biom <- tryCatch(as_rbiom(biom), error = function (e) NULL)

  if (!inherits(biom, 'rbiom'))
    cli_abort("Can't apply metadata partitions to a non-rbiom object.")

  return (biom)
}



#________________________________________________________
# The subset of `biom` for one stratum of `vars`, with a history entry
# such as `subset(biom, Sex == "Male" & `Body Site` == "Stool")`.
#________________________________________________________
ply_subset <- function (biom, df, vars) {

  conds <- lapply(vars, function (v) {
    val <- df[[v]][[1]]
    if (is.factor(val)) val <- as.character(val)
    call('==', as.name(v), val)
  })

  pc <- prov_call('subset', list(quote(biom), Reduce(function (a, b) call('&', a, b), conds)))

  sub_biom <- biom$clone()

  biom_transform(
    biom         = sub_biom,
    expr         = sub_biom$counts <- sub_biom$counts[, df[['.sample']], drop = FALSE],
    call         = pc$call,
    reproducible = pc$reproducible,
    clone        = TRUE )
}
