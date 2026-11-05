#' Print code blocks neatly.
#'
#' @param x   An object of class `rbiom_code`.
#' @param ...   Not used.
#'
#' @noRd
#' @keywords internal
#' @export
print.rbiom_code <- function (x, ...) {

  x <- paste(unclass(x), collapse = "\n")

  if (nzchar(system.file(package = "prettycode")))
    x <- strsplit(x, '\n')[[1]] %>%
      prettycode::highlight() %>%
      paste(collapse = '\n')

  cat(x, "\n", sep = "")

  return (invisible(NULL))
}



#' Also search attributes.
#'
#' `$code` assembles the full reproducible script (library calls, the rbiom
#' object's history, the data step, and the statistics / ggplot2 code).
#'
#' @noRd
#' @keywords internal
#' @export
`$.rbiom_tbl` <- function (obj, nm) {

  if (!nm %in% c('cmd', 'code', 'history', 'stats', 'taxa_coords', 'taxa_stats')) {
    NextMethod()

  } else if (hasName(obj, nm)) {
    NextMethod()

  } else if (nm == 'code') {
    assemble_code(obj)

  } else {

    val <- attr(obj, nm, exact = TRUE)

    if (is.null(val) && hasName(obj, 'data'))
      val <- attr(obj[['data']], nm, exact = TRUE)

    if (!is.null(val) && nm %in% c('cmd', 'history'))
      return (add_class(paste(val, collapse = "\n"), 'rbiom_code'))

    return (val)
  }

}


#' @export
`$.rbiom_plot` <- `$.rbiom_tbl`



#' Adds attr(,'tbl_sum') to tibble header.
#'
#' @noRd
#' @keywords internal
#' @export
tbl_sum.rbiom_tbl <- function (x) {
  c(attr(x, 'tbl_sum'), NextMethod())
}





#' Convert an rbiom object to a base R list.
#'
#' @inherit documentation_biom.rbiom
#'
#' @family conversion
#'
#' @param ...   Not used.
#'
#' @return A list with names
#'         `c('counts', 'metadata', 'taxonomy', 'tree', 'sequences', 'id', 'comment', 'date', 'generated_by')`.
#'
#' @export
#'
as.list.rbiom <- function (x, ...) {
  biom <- x
  list(
    'counts'       = biom$counts,
    'metadata'     = biom$metadata,
    'taxonomy'     = biom$taxonomy,
    'tree'         = biom$tree,
    'sequences'    = biom$sequences,
    'id'           = biom$id,
    'comment'      = biom$comment,
    'date'         = biom$date,
    'generated_by' = biom$generated_by )
}





#' Convert an rbiom object to a simple count matrix.
#'
#' Identical to running `as.matrix(biom$counts)`.
#'
#' @inherit documentation_biom.rbiom
#'
#' @family conversion
#'
#' @param ...   Not used.
#'
#' @return A base R matrix with OTUs as rows and samples as columns.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     as.matrix(hmp50)[1:5,1:5]
#'
as.matrix.rbiom <- function (x, ...) {
  as.matrix(x$counts)
}




#' Map sample names to metadata field values.
#'
#' @inherit documentation_biom.rbiom
#'
#' @family samples
#'
#' @param var   The metadata field name specified as:
#' \itemize{
#'   \item{The metadata field name to retrieve. Can be abbreviated.}
#'   \item{A positive integer, giving the position counting from the left.}
#'   \item{A negative integer, giving the position counting from the right.}
#' }
#' Default: `-1`
#'
#' @param name   The column to be used as names for a named vector.
#'        Specified in a similar manner as var. Default: `".sample"`
#'
#' @param ...   Not used.
#'
#' @return A vector of metadata values, named with sample names.
#'
#' @seealso `taxa_map()`
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     pull(hmp50, 'Age') %>% head()
#'
#'     pull(hmp50, 'bod') %>% head(4)
#'
pull.rbiom <- function (.data, var = -1, name = ".sample", ...) {
  biom <- .data
  if (!is_integerish(var))  validate_biom_field("var")
  if (!is_integerish(name)) validate_biom_field("name", null_ok = TRUE)
  dplyr::pull(.data = biom$metadata, var = var, name = name)
}




#' Get a glimpse of your metadata.
#'
#' @inherit documentation_biom.rbiom
#'
#' @family metadata
#'
#' @param width   Width of output. See [tibble::glimpse()] documentation.
#'        Default: `NULL`
#'
#' @param ...   Not used.
#'
#' @return The original `biom`, invisibly.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     glimpse(hmp50)
#'
glimpse.rbiom <- function (x, width = NULL, ...) {
  biom <- x
  eval.parent(tibble::glimpse(x = biom$metadata, width = width))
  return (invisible(biom))
}




#________________________________________________________
# Shared by the transformation methods below: run `expr` (which edits
# `biom` through its setters) without the setters logging, then log the
# user-level call `call` instead.
#________________________________________________________
biom_transform <- function (biom, expr, call, reproducible = TRUE, note = NULL, clone = TRUE) {
  prov_quiet(expr)
  prov_log(biom, call = call, reproducible = reproducible, note = note)
  if (isTRUE(clone)) { return (biom) } else { return (invisible(biom)) }
}




#' Create, modify, and delete metadata fields.
#'
#' mutate() creates new fields in `$metadata` that are functions of existing
#' metadata fields. It can also modify (if the name is the same as an existing
#' field) and delete fields (by setting their value to NULL).
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#'
#' @name modify_metadata
#' @family transformations
#'
#' @param ...   Passed on to [dplyr::mutate()] or [dplyr::rename()].
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- slice_max(hmp50, BMI, n = 6)
#'     biom$metadata
#'
#'     # Add a new field to the metadata
#'     biom <- mutate(biom, Obsese = BMI >= 30)
#'     biom$metadata
#'
#'     # Rename a metadata field
#'     biom <- rename(biom, 'Age (years)' = "Age")
#'     biom$metadata
#'
mutate.rbiom <- function (.data, ..., clone = TRUE) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  env  <- parent.frame()

  pc <- prov_inline_dots(rlang::enexprs(...), env, exclude = biom$fields)

  biom_transform(
    biom         = biom,
    expr         = biom$metadata <- eval.parent(dplyr::mutate(.data = biom$metadata, ...)),
    call         = as.call(c(quote(mutate), quote(biom), pc$exprs)),
    reproducible = pc$reproducible,
    clone        = clone )
}



#' @rdname modify_metadata
#' @export
rename.rbiom <- function (.data, ..., clone = TRUE) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  env  <- parent.frame()

  pc <- prov_inline_dots(rlang::enexprs(...), env, exclude = biom$fields)

  biom_transform(
    biom         = biom,
    expr         = biom$metadata <- eval.parent(dplyr::rename(.data = biom$metadata, ...)),
    call         = as.call(c(quote(rename), quote(biom), pc$exprs)),
    reproducible = pc$reproducible,
    clone        = clone )
}


#________________________________________________________
# Inline free variables in each of a list of captured expressions.
#________________________________________________________
prov_inline_dots <- function (exprs, env, exclude = character()) {

  ok <- TRUE

  exprs <- lapply(exprs, function (e) {
    if (!is.language(e)) return (e)
    pi <- prov_inline(e, env, exclude)
    if (!pi$reproducible) ok <<- FALSE
    pi$expr
  })

  list(exprs = exprs, reproducible = ok)
}





#' Evaluate expressions on metadata.
#'
#' `with()` will return the result of your expression. `within()` will return
#' an rbiom object.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_default
#'
#' @name with
#' @family transformations
#'
#' @param expr   Passed on to [base::with()] or [base::within()].
#'
#' @param ...   Not used.
#'
#' @return See description.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     with(hmp50, table(`Body Site`, Sex))
#'
#'     biom <- within(hmp50, {
#'       age_bin = cut(Age, 5)
#'       bmi_bin = cut(BMI, 5)
#'     })
#'     biom$metadata
#'
with.rbiom <- function (data, expr, ...) {
  eval(expr = substitute(expr), envir = data$metadata, enclos = parent.frame())
}



#' @rdname with
#' @export
within.rbiom <- function (data, expr, clone = TRUE, ...) {

  biom   <- if (isTRUE(clone)) data$clone() else data
  parent <- parent.frame()
  expr   <- substitute(expr)
  pi     <- prov_inline(expr, parent, exclude = biom$fields)

  biom_transform(
    biom = biom,
    expr = {

      data <- biom$metadata

      # copied from base::within.data.frame
      e <- evalq(environment(), data, parent)
      eval(expr, e)
      l <- as.list(e, all.names = TRUE)
      l <- l[!vapply(l, is.null, NA, USE.NAMES = FALSE)]
      nl <- names(l)
      del <- setdiff(names(data), nl)
      data[nl] <- l
      data[del] <- NULL

      biom$metadata <- data
    },
    call         = call('within', quote(biom), pi$expr),
    reproducible = pi$reproducible,
    clone        = clone )
}





#' Subset an rbiom object by sample names, OTU names, metadata, or taxonomy.
#'
#' Dropping samples or OTUs will lead to observations being removed from the
#' OTU matrix (`biom$counts`). OTUs and samples with zero observations are
#' automatically removed from the rbiom object.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#'
#' @name subset
#' @family transformations
#'
#' @param subset   Logical expression for rows to keep. See [base::subset()].
#'        As with [base::subset()], `NA` values are treated as `FALSE`.
#'
#' @param i,j   The sample or OTU names to keep. Or a logical/integer vector
#'        indicating which sample names from `biom$samples` or `biom$otus` to
#'        keep. Subsetting with `[i]` takes `i` as samples, whereas `[i,j]`
#'        takes `i` as otus and `j` as samples (corresponding to `[rows, cols]`
#'        in the underlying `biom$counts` matrix).
#'
#' @param fields   Which metadata field(s) to check for `NA`s, or `".all"` to
#'        check all metadata fields.
#'
#' @param drop   Not used
#'
#' @param ...   Not used.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'     library(dplyr)
#'
#'     # Subset to specific samples
#'     biom <- hmp50[c('HMP20', 'HMP42', 'HMP12')]
#'     biom$metadata
#'
#'     # Subset to specific OTUs
#'     biom <- hmp50[c('LtbAci52', 'UncO2012'),] # <- Trailing ,
#'     biom$taxonomy
#'
#'     # Subset to specific samples and OTUs
#'     biom <- hmp50[c('LtbAci52', 'UncO2012'), c('HMP20', 'HMP42', 'HMP12')]
#'     as.matrix(biom)
#'
#'     # Subset samples according to metadata
#'     biom <- subset(hmp50, `Body Site` %in% c('Saliva') & Age < 25)
#'     biom$metadata
#'
#'     # Subset OTUs according to taxonomy
#'     biom <- subset_taxa(hmp50, Phylum == 'Cyanobacteria')
#'     biom$taxonomy
#'
#'     # Remove samples with NA metadata values
#'     biom <- mutate(hmp50, BS2 = na_if(`Body Site`, 'Saliva'))
#'     biom$metadata
#'     biom <- na.omit(biom)
#'     biom$metadata
#'
#'     # The steps that produced an object are recorded.
#'     biom$code
#'

subset.rbiom <- function (x, subset, clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) x$clone() else x
  expr <- substitute(subset)
  env  <- parent.frame()

  keep <- eval(expr = expr, envir = biom$metadata, enclos = env)

  if (is.logical(keep) && length(keep) == 1)
    keep <- rep(keep, biom$n_samples)

  if (!is.logical(keep) || length(keep) != biom$n_samples)
    cli_abort(c(
      'x' = "`subset` must evaluate to a logical vector with one value per sample.",
      'i' = "Got {.type {keep}} of length {length(keep)}; expected length {biom$n_samples}." ))

  keep[is.na(keep)] <- FALSE

  if (!any(keep))
    cli_abort("No samples match `subset`.")

  pi <- prov_inline(expr, env, exclude = biom$fields)

  biom_transform(
    biom         = biom,
    expr         = biom$counts <- biom$counts[, keep, drop = FALSE],
    call         = call('subset', quote(biom), pi$expr),
    reproducible = pi$reproducible,
    clone        = clone )
}


#' @rdname subset
#' @export
`[.rbiom` <- function (x, i, j, ..., clone = TRUE, drop = FALSE) {

  biom <- if (isTRUE(clone)) x$clone() else x

  # NA in a logical index means "drop", as in base::subset().
  na_false <- function (idx, what, choices) {
    if (is.logical(idx)) idx[is.na(idx)] <- FALSE
    if (is.character(idx) && length(bad <- setdiff(idx, choices)) > 0)
      cli_abort("Unknown {what} name{?s}: {.val {bad}}.")
    idx
  }

  # Number of index positions in the call: biom[i] vs biom[i, j].
  n_idx <- nargs() - 1L - (!missing(clone)) - (!missing(drop)) - length(list(...))

  if (n_idx < 2) {
    # Called as biom[samples]
    j <- na_false(i, "sample", biom$samples)
    pc <- prov_call('[', list(quote(biom), j))
    biom_transform(
      biom         = biom,
      expr         = biom$counts <- biom$counts[, j, drop = FALSE],
      call         = pc$call,
      reproducible = pc$reproducible,
      clone        = clone )

  } else {
    # Called as biom[otus, samples], biom[otus, ], or biom[, samples]
    idx <- alist(, )
    if (!missing(i)) { i <- na_false(i, "OTU",    biom$otus);    idx[1] <- list(i) }
    if (!missing(j)) { j <- na_false(j, "sample", biom$samples); idx[2] <- list(j) }
    pc <- prov_call('[', c(list(quote(biom)), idx))

    biom_transform(
      biom = biom,
      expr = {
        if      (missing(i)) { biom$counts <- biom$counts[ , j, drop = FALSE] }
        else if (missing(j)) { biom$counts <- biom$counts[i,  , drop = FALSE] }
        else                 { biom$counts <- biom$counts[i, j, drop = FALSE] }
      },
      call         = pc$call,
      reproducible = pc$reproducible,
      clone        = clone )
  }
}




#' @rdname subset
#' @export
na.omit.rbiom <- function (object, fields = ".all", clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) object$clone() else object

  if (length(fields) == 0) return (biom)

  if (eq(fields, ".all")) fields <- biom$fields
  validate_biom_field('fields', max = Inf)

  keep <- stats::complete.cases(biom$metadata[,fields,drop=FALSE])

  if (!any(keep))
    cli_abort("All samples have NA values in {.field {fields}}.")

  args <- if (identical(fields, biom$fields)) list() else list(fields = fields)

  biom_transform(
    biom  = biom,
    expr  = biom$counts <- biom$counts[, keep, drop = FALSE],
    call  = as.call(c(quote(na.omit), quote(biom), args)),
    clone = clone )
}





#' Subset to a specific number of samples.
#'
#' @inherit documentation_biom.rbiom
#' @inherit documentation_return.biom return
#' @inherit documentation_default
#' @inheritParams dplyr::slice
#'
#' @name slice_metadata
#' @family transformations
#'
#' @param ...   For `slice()`, integer row indexes. For other `slice_*()`
#'        functions, not used. See [dplyr::slice()].
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     # The last 3 samples in the metadata table.
#'     biom <- slice_tail(hmp50, n = 3)
#'     biom$metadata
#'
#'     # The 3 oldest subjects sampled.
#'     biom <- slice_max(hmp50, Age, n = 3)
#'     biom$metadata
#'
#'     # Pick 3 samples at random.
#'     biom <- slice_sample(hmp50, n = 3)
#'     biom$metadata
#'
slice.rbiom <- function (.data, ..., .by = NULL, .preserve = FALSE, clone = TRUE) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  env  <- parent.frame()

  pc   <- prov_inline_dots(rlang::enexprs(...), env, exclude = biom$fields)
  args <- pc$exprs
  if (!is.null(.by_expr <- substitute(.by))) args[['.by']] <- .by_expr
  if (!isFALSE(.preserve))                   args[['.preserve']] <- .preserve

  biom_transform(
    biom = biom,
    expr = {
      df <- biom$metadata
      df <- eval.parent(dplyr::slice(df, ..., .by = {{.by}}, .preserve = .preserve))
      suppressWarnings(biom$metadata <- df)
    },
    call         = as.call(c(quote(slice), quote(biom), args)),
    reproducible = pc$reproducible,
    clone        = clone )
}


#________________________________________________________
# Collect the non-missing, non-default arguments of a slice_*() call.
#________________________________________________________
slice_args <- function (n, prop, by_expr, ...) {
  args <- list()
  if (!missing(n))        args[['n']]    <- n
  if (!missing(prop))     args[['prop']] <- prop
  if (!is.null(by_expr))  args[['by']]   <- by_expr
  c(args, list(...))
}


#' @rdname slice_metadata
#' @export
slice_head.rbiom <- function (.data, n, prop, by = NULL, clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  args <- slice_args(n, prop, substitute(by))

  biom_transform(
    biom = biom,
    expr = {
      df <- biom$metadata
      df <- eval.parent(dplyr::slice_head(df, n = n, prop = prop, by = {{by}}))
      suppressWarnings(biom$metadata <- df)
    },
    call  = as.call(c(quote(slice_head), quote(biom), args)),
    clone = clone )
}


#' @rdname slice_metadata
#' @export
slice_tail.rbiom <- function (.data, n, prop, by = NULL, clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  args <- slice_args(n, prop, substitute(by))

  biom_transform(
    biom = biom,
    expr = {
      df <- biom$metadata
      df <- eval.parent(dplyr::slice_tail(df, n = n, prop = prop, by = {{by}}))
      suppressWarnings(biom$metadata <- df)
    },
    call  = as.call(c(quote(slice_tail), quote(biom), args)),
    clone = clone )
}


#' @rdname slice_metadata
#' @export
slice_min.rbiom <- function (.data, order_by, n, prop, by = NULL, with_ties = TRUE, na_rm = FALSE, clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  args <- c(
    list(order_by = substitute(order_by)),
    slice_args(n, prop, substitute(by)),
    if (!isTRUE(with_ties)) list(with_ties = with_ties),
    if (!isFALSE(na_rm))    list(na_rm = na_rm) )

  biom_transform(
    biom = biom,
    expr = {
      df <- eval.parent(
        expr = dplyr::slice_min(
          .data     = biom$metadata,
          order_by  = {{order_by}},
          n         = n,
          prop      = prop,
          by        = {{by}},
          with_ties = with_ties,
          na_rm     = na_rm ))
      suppressWarnings(biom$metadata <- df)
    },
    call  = as.call(c(quote(slice_min), quote(biom), args)),
    clone = clone )
}


#' @rdname slice_metadata
#' @export
slice_max.rbiom <- function (.data, order_by, n, prop, by = NULL, with_ties = TRUE, na_rm = FALSE, clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  args <- c(
    list(order_by = substitute(order_by)),
    slice_args(n, prop, substitute(by)),
    if (!isTRUE(with_ties)) list(with_ties = with_ties),
    if (!isFALSE(na_rm))    list(na_rm = na_rm) )

  biom_transform(
    biom = biom,
    expr = {
      df <- eval.parent(
        expr = dplyr::slice_max(
          .data     = biom$metadata,
          order_by  = {{order_by}},
          n         = n,
          prop      = prop,
          by        = {{by}},
          with_ties = with_ties,
          na_rm     = na_rm ))
      suppressWarnings(biom$metadata <- df)
    },
    call  = as.call(c(quote(slice_max), quote(biom), args)),
    clone = clone )
}


#' @rdname slice_metadata
#' @export
slice_sample.rbiom <- function (.data, n, prop, by = NULL, weight_by = NULL, replace = FALSE, clone = TRUE, ...) {

  biom <- if (isTRUE(clone)) .data$clone() else .data
  args <- c(
    slice_args(n, prop, substitute(by)),
    if (!is.null(wb <- substitute(weight_by))) list(weight_by = wb),
    if (!isFALSE(replace)) list(replace = replace) )

  biom_transform(
    biom = biom,
    expr = {
      df <- biom$metadata
      df <- eval.parent(dplyr::slice_sample(df, n = n, prop = prop, by = {{by}},  weight_by = {{weight_by}}, replace = replace))
      suppressWarnings(biom$metadata <- df)
    },
    call  = as.call(c(quote(slice_sample), quote(biom), args)),
    note  = "random; depends on the RNG state",
    clone = clone )
}
