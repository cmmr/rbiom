#________________________________________________________
# Split-apply-combine helpers. These replace the subset of plyr that rbiom
# used (ddply, dlply, adply, ldply, empty) with base R and dplyr::bind_rows().
#________________________________________________________


#________________________________________________________
# Coerce whatever a user function returned into data.frame rows, or NULL
# when it returned nothing. Named atomic vectors become a single row.
#________________________________________________________
as_rows <- function (x) {

  if (is.null(x))       return (NULL)
  if (is.data.frame(x)) return (if (nrow(x) == 0) NULL else x)
  if (is.matrix(x))     return (if (nrow(x) == 0) NULL else as.data.frame(x, optional = TRUE))

  if (is.atomic(x)) {
    if (length(x) == 0) return (NULL)
    if (is.null(names(x))) names(x) <- paste0('V', seq_along(x))
    return (as.data.frame(as.list(x), optional = TRUE, stringsAsFactors = FALSE))
  }

  if (is.list(x)) {
    if (length(x) == 0) return (NULL)
    if (is.null(names(x))) names(x) <- paste0('V', seq_along(x))
    return (as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE))
  }

  cli_abort("Can't convert {.type {x}} to data.frame rows.")
}


#________________________________________________________
# Apply `fun` to each subset of `df` defined by the unique combinations of
# the `by` columns.
#
#   as = "df":   row-bind the results, prepending the `by` values
#                (like plyr::ddply).
#   as = "list": return a list of results with a 'split_labels' data.frame
#                attribute (like plyr::dlply).
#
# Groups are ordered by sorted `by` values (factor level order). Rows with
# an NA in any `by` column are dropped. Empty combinations are dropped.
#________________________________________________________
split_apply <- function (df, by = NULL, fun, ..., as = "df") {

  stopifnot(is.data.frame(df), is.function(fun))
  as <- match.arg(as, c("df", "list"))

  if (length(by) == 0) {
    res <- fun(df, ...)
    if (as == "list") return (structure(list(res), split_labels = data.frame()))
    return (as_tibble(as_rows(res) %||% data.frame()))
  }

  if (length(missing <- setdiff(by, names(df))) > 0)
    cli_abort("Can't split by missing column{?s}: {.val {missing}}.")

  keys <- df[, by, drop = FALSE]
  f    <- lapply(keys, function (k) if (is.factor(k)) droplevels(k) else factor(k))
  idx  <- split(seq_len(nrow(df)), f, drop = TRUE, lex.order = TRUE)

  labels <- keys[vapply(idx, `[[`, 1L, 1L), , drop = FALSE]
  labels <- as.data.frame(labels, optional = TRUE)
  rownames(labels) <- NULL

  results <- lapply(idx, function (i) fun(df[i, , drop = FALSE], ...))

  if (as == "list") {
    names(results) <- names(idx)
    attr(results, 'split_labels') <- labels
    return (results)
  }

  rows <- lapply(seq_along(results), function (j) {
    r <- as_rows(results[[j]])
    if (is.null(r)) return (NULL)
    cbind(labels[rep(j, nrow(r)), , drop = FALSE], r)
  })

  res <- dplyr::bind_rows(rows)
  if (nrow(res) == 0) res <- labels[0, , drop = FALSE]

  return (as_tibble(res))
}


#________________________________________________________
# Apply `fun` to each row of `df` (passed as a one-row data.frame) and
# row-bind the results, prepending that row's values (like plyr::adply).
#________________________________________________________
rows_apply <- function (df, fun, ...) {

  stopifnot(is.data.frame(df), is.function(fun))

  rows <- lapply(seq_len(nrow(df)), function (i) {
    row <- df[i, , drop = FALSE]
    r   <- as_rows(fun(row, ...))
    if (is.null(r)) return (NULL)
    rownames(row) <- NULL
    cbind(row[rep(1L, nrow(r)), , drop = FALSE], r)
  })

  res <- dplyr::bind_rows(rows)
  if (nrow(res) == 0) res <- df[0, , drop = FALSE]

  return (as_tibble(res))
}


#________________________________________________________
# Is a data.frame/matrix/vector empty? (Replaces plyr::empty.)
#________________________________________________________
is_empty_df <- function (x) {
  is.null(x) || length(x) == 0 || NROW(x) == 0 || (is.data.frame(x) && ncol(x) == 0)
}
