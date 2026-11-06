#________________________________________________________
# Shared helpers for the stats_* functions.
#
# Each statistics function has two parallel parts that are written side by
# side so they stay in step:
#
#   * closures that compute the result from the real, evaluated arguments
#     (stats_run() below provides the NA handling, per-stratum and
#     per-pair scaffolding, and the p-value adjustment), and
#   * the R code text returned to the user as `$code`, which records the
#     same computation. That text is a record only: the package never
#     parses or evaluates it.
#
# Anything that originates from the data (column names, group labels) is
# quoted when it is written into the code text: names with coan() and
# values with double_quote() / deparse_vec(), so that hostile names can't
# change the meaning of the code.
#
# The generated code is written to run with only `library(rbiom)` attached
# (rbiom re-exports `%>%` and the dplyr verbs it uses); everything else is
# referenced with `pkg::`.
#________________________________________________________



#________________________________________________________
# Validate the stats arguments shared by stats_table() and friends, and
# report how many rows will be dropped because of NA values. Returns the
# complete-case subset of the relevant columns (invisibly).
#________________________________________________________
stats_validate <- function (df, regr, resp, stat.by, split.by, strata) {

  cols <- unique(c(regr, resp, stat.by, split.by, strata))

  if (as.logical(anyDuplicated(c(regr, resp, stat.by, split.by, strata))))
    cli_abort("`regr`, `resp`, `stat.by`, `split.by`, and `strata` must all be unique.")

  sub  <- df[, cols, drop = FALSE]
  keep <- stats::complete.cases(sub)

  if (!any(keep))
    cli_abort("No rows remain after removing rows with NA values in {.field {cols}}.")

  if (!all(keep))
    cli::cli_inform(c(
      'i' = "Dropped {sum(!keep)} of {length(keep)} row{?s} with NA values in {.field {cols}}." ))

  sub <- sub[keep, , drop = FALSE]


  #________________________________________________________
  # Group sizes: pairwise/group-wise tests need two or more groups in
  # every stratum, and each group needs at least two observations.
  #________________________________________________________
  if (!is.null(stat.by)) {

    groups <- as.character(sub[[stat.by]])

    if (length(unique(groups)) < 2)
      cli_abort(c(
        'x' = "`stat.by` field {.field {stat.by}} has only one group ({.val {unique(groups)}}).",
        'i' = "At least two groups are needed to compare." ))

    # A stratum with a single group gets an NA row and is named in the
    # warning that stats_report() raises.
    strata_id <- if (length(split.by) == 0) rep("", nrow(sub)) else
      do.call(paste, c(lapply(sub[, split.by, drop = FALSE], as.character), sep = " / "))

    sizes <- table(paste(strata_id, groups, sep = " / "))
    small <- sub("^ / ", "", names(sizes)[sizes < 2])
    if (length(small) > 0)
      cli_warn(c(
        '!' = "{qty(length(small))}Group{?s} with fewer than two observations: {.val {small}}.",
        'i' = "Tests involving {qty(length(small))}{?this/these} group{?s} will fail or be unreliable." ))
  }

  return (invisible(sub))
}



#________________________________________________________
# Run a test on `df`.
#
#   pre  - NULL, or list(fun = function (data) -> value, code = text):
#          computed once per stratum and passed to `test$fun` as `at`.
#   test - list(fun = function (data, pair, at) -> tibble, code = text).
#          For pairwise tests, `data` is the pair's subset and `pair` the
#          two group labels; otherwise `pair` is NULL.
#
# The structure here mirrors stats_wrap() line for line: the same column
# selection, NA filtering, factor handling, per-pair subsetting, error
# handling, and p-value adjustment.
#________________________________________________________
stats_run <- function (df, regr, resp, stat.by, split.by, strata, pre, test, pairwise, p.adj) {

  #________________________________________________________
  # Prologue: standardized column names, complete cases.
  #________________________________________________________
  cols <- c(
    .resp = resp,
    if (!is.null(regr))    c(.regr    = regr),
    if (!is.null(stat.by)) c(.stat.by = stat.by),
    if (!is.null(strata))  c(.strata  = strata),
    stats::setNames(split.by, split.by) )

  data <- dplyr::select(df, !!!rlang::syms(cols))
  data <- data[stats::complete.cases(data), , drop = FALSE]
  if (!is.null(stat.by)) data$.stat.by <- droplevels(as.factor(data$.stat.by))
  if (!is.null(strata))  data$.strata  <- droplevels(as.factor(data$.strata))
  for (col in split.by)  data[[col]]   <- droplevels(as.factor(data[[col]]))

  multi <- !is.null(split.by) || isTRUE(pairwise)

  na_row <- function (pair) {
    if (is.null(pair)) return (tibble::tibble(.p.val = NA_real_))
    tibble::tibble(!!stat.by := paste(pair, collapse = ' - '), .p.val = NA_real_)
  }


  #________________________________________________________
  # One stratum: every pair of groups, or a single test.
  #________________________________________________________
  stratum <- function (data) {

    at <- if (!is.null(pre)) pre$fun(data)

    if (isTRUE(pairwise)) {

      lvls <- levels(droplevels(data$.stat.by))
      if (length(lvls) < 2) stop("fewer than two groups in this stratum")
      pairs <- combn(lvls, 2, simplify = FALSE)
      dplyr::bind_rows(lapply(pairs, function (pair) {
        tryCatch(
          error = function (e) { warning(conditionMessage(e), call. = FALSE); na_row(pair) },
          expr  = {
            pdata <- data[data$.stat.by %in% pair, , drop = FALSE]
            test$fun(pdata, pair, at)
          })
      }))

    } else {
      test$fun(data, NULL, at)
    }
  }


  #________________________________________________________
  # Split.by strata (each guarded, so one failing stratum becomes an NA
  # row), or a single evaluation.
  #________________________________________________________
  stats <- if (!is.null(split.by)) {
    data %>%
      dplyr::group_by(!!!rlang::syms(split.by)) %>%
      dplyr::group_modify(function (data, keys) {
        tryCatch(
          error = function (e) { warning(conditionMessage(e), call. = FALSE); na_row(NULL) },
          expr  = stratum(data) )
      }) %>%
      dplyr::ungroup()
  } else {
    stratum(data)
  }


  #________________________________________________________
  # Correct for multiple comparisons and arrange by p.val.
  #________________________________________________________
  .p.val <- NULL # for CRAN check only
  stats <- dplyr::mutate(stats, .adj.p = p.adjust(.p.val, p.adj), .after = .p.val)
  stats <- dplyr::arrange(stats, .p.val)

  return (stats)
}



#________________________________________________________
# Code formatting helpers.
#________________________________________________________

# Indent all but the first line of `code` by `n` spaces.
indent_code <- function (code, n = 2) {
  gsub("\n", paste0("\n", strrep(" ", n)), code, fixed = TRUE)
}

# Format `tibble::tibble(name = expr, ...)` with aligned (quoted) names.
tibble_code <- function (cols, fn = "tibble::tibble") {
  cols <- unlist(cols)
  nms  <- coan(names(cols))
  fmt  <- paste0("%-", max(nchar(nms)), "s = %s")
  paste0(fn, "(\n  ", paste(sprintf(fmt, nms, cols), collapse = ",\n  "), " )")
}


#________________________________________________________
# A literal R expression that recreates an atomic vector (names, factor
# levels, and full double precision are preserved), for embedding data in
# generated code.
#________________________________________________________
deparse_vec <- function (x, indent = 2) {

  # Doubles are written so that they parse back exactly on every platform:
  # see num_literal() in provenance.r.
  opts <- c("keepNA", "keepInteger", "niceNames")
  sep  <- paste0("\n", strrep(" ", indent))
  dp   <- function (v) paste(deparse_exact(v, width.cutoff = 100L, control = opts), collapse = sep)

  if (is.factor(x)) {
    vals <- stats::setNames(as.character(x), names(x))
    return (sprintf("factor(%s,%slevels = %s)", dp(vals), sep, dp(levels(x))))
  }

  dp(stats::setNames(as.vector(x), names(x)))
}


#________________________________________________________
# Attach provenance to a result that is not a data.frame (a matrix, a
# dist, a vector) and give it the 'rbiom_obj' class so that `$cmd`,
# `$history`, and `$code` work on it as they do on tibbles. S4 objects
# (sparse matrices) only get the attributes.
#________________________________________________________
prov_result <- function (res, cmd = NULL, biom = NULL, data = NULL) {

  res <- prov_attach(res, cmd, biom, data)

  if (!is.null(res) && !is.data.frame(res) && !isS4(res) && !inherits(res, 'rbiom_obj'))
    class(res) <- c('rbiom_obj', class(res))

  return (res)
}



#' Provenance accessors for non-data.frame results.
#'
#' `$cmd`, `$history`, and `$code` on the matrices, distance matrices, and
#' vectors returned by rbiom functions. Printing hides the provenance
#' attributes.
#'
#' @noRd
#' @keywords internal
#' @export
`$.rbiom_obj` <- function (obj, nm) {

  if (nm == 'code') return (assemble_code(obj))

  val <- attr(obj, nm, exact = TRUE)

  if (!is.null(val) && nm %in% c('cmd', 'history'))
    return (add_class(paste(val, collapse = "\n"), 'rbiom_code'))

  return (val)
}


#' @export
print.rbiom_obj <- function (x, ...) {

  y <- x
  for (a in c('cmd', 'code', 'history', 'data_cmd'))
    attr(y, a) <- NULL

  cls      <- setdiff(oldClass(y), 'rbiom_obj')
  implicit <- class(unclass(y))
  oldClass(y) <- if (length(cls) > 0 && !identical(cls, implicit)) cls else NULL

  print(y, ...)

  return (invisible(x))
}



#________________________________________________________
# The name of a p.adjust() method for captions, or NULL for "none".
#________________________________________________________
p_adj_label <- function (method) {
  switch(
    EXPR       = method,
    holm       = "Holm",
    hochberg   = "Hochberg",
    hommel     = "Hommel",
    BH         = "FDR (Benjamini-Hochberg)",
    fdr        = "FDR (Benjamini-Hochberg)",
    BY         = "FDR (Benjamini-Yekutieli)",
    bonferroni = "Bonferroni",
    none       = NULL,
    method )
}



#________________________________________________________
# A copy of a function's parameters for current_cmd(), in which objects
# that can't be deparsed (a phylo tree) are shown by name.
#________________________________________________________
cmd_params <- function (params) {
  p <- as.list(params)
  if (inherits(p$tree, 'phylo')) p$tree <- aa(NA, display = "tree")
  return (p)
}


#________________________________________________________
# Remove provenance attributes, e.g. before hashing a data.frame for the
# cache so that identical data from different pipelines share an entry.
#________________________________________________________
strip_prov <- function (x) {
  for (a in c('cmd', 'code', 'history', 'data_cmd'))
    attr(x, a) <- NULL
  return (x)
}



#________________________________________________________
# The code that selects, renames, and NA-filters the stats columns
# (mirrors the prologue of stats_run()).
#________________________________________________________
stats_prologue <- function (regr, resp, stat.by, split.by, strata) {

  sel <- c(
    sprintf(".resp = %s", coan(resp)),
    if (!is.null(regr))    sprintf(".regr = %s",    coan(regr)),
    if (!is.null(stat.by)) sprintf(".stat.by = %s", coan(stat.by)),
    if (!is.null(strata))  sprintf(".strata = %s",  coan(strata)),
    if (!is.null(split.by)) coan(split.by) )

  c(
    sprintf("data <- dplyr::select(data, %s)", paste(sel, collapse = ", ")),
    "data <- data[stats::complete.cases(data), , drop = FALSE]",
    if (!is.null(stat.by)) "data$.stat.by <- droplevels(as.factor(data$.stat.by))",
    if (!is.null(strata))  "data$.strata <- droplevels(as.factor(data$.strata))",
    vapply(split.by, FUN.VALUE = '', function (col) {
      sprintf("data[[%s]] <- droplevels(as.factor(data[[%s]]))", double_quote(col), double_quote(col)) }) ) %>%
    paste(collapse = "\n")
}



#________________________________________________________
# Wrap the per-stratum test code in the pairwise / split.by scaffolding
# (mirrors stratum() in stats_run()).
#________________________________________________________
stats_wrap <- function (pre_code, test_code, stat.by, split.by, pairwise) {


  #________________________________________________________
  # Per-pair function (each pair guarded), or a single test.
  #________________________________________________________
  if (isTRUE(pairwise)) {

    na_row <- sprintf(
      "tibble::tibble(%s = paste(pair, collapse = ' - '), .p.val = NA_real_)",
      coan(stat.by) )

    inner <- paste(sep = "\n",
      "lvls  <- levels(droplevels(data$.stat.by))",
      "if (length(lvls) < 2) stop('fewer than two groups in this stratum')",
      "pairs <- combn(lvls, 2, simplify = FALSE)",
      "dplyr::bind_rows(lapply(pairs, function (pair) {",
      "  tryCatch(",
      sprintf("    error = function (e) { warning(conditionMessage(e), call. = FALSE); %s },", na_row),
      "    expr  = {",
      "      pdata <- data[data$.stat.by %in% pair, , drop = FALSE]",
      paste0("      ", indent_code(test_code, 6)),
      "    })",
      "}))" )

  } else {
    inner <- test_code
  }

  if (!is.null(pre_code))
    inner <- paste(pre_code, "", inner, sep = "\n")


  #________________________________________________________
  # Split.by strata (each guarded), or a single evaluation.
  #________________________________________________________
  if (!is.null(split.by)) {

    paste(sep = "\n",
      "stats <- data %>%",
      sprintf("  dplyr::group_by(%s) %%>%%", paste(coan(split.by), collapse = ", ")),
      "  dplyr::group_modify(function (data, keys) {",
      "    tryCatch(",
      "      error = function (e) { warning(conditionMessage(e), call. = FALSE); tibble::tibble(.p.val = NA_real_) },",
      "      expr  = {",
      paste0("        ", indent_code(inner, 8)),
      "      })",
      "  }) %>%",
      "  dplyr::ungroup()" )

  } else {

    paste(sep = "\n",
      "stats <- local({",
      paste0("  ", indent_code(inner, 2)),
      "})" )
  }
}



#________________________________________________________
# Run the test, report failures in one message, and decorate the result
# with its code and header. `model` is the text shown as the tibble header.
#________________________________________________________
stats_finalize <- function (df, pre, test, regr, resp, stat.by, split.by, strata, p.adj, pairwise, model) {

  code <- paste(sep = "\n\n",
    stats_prologue(regr, resp, stat.by, split.by, strata),
    stats_wrap(pre$code, test$code, stat.by, split.by, pairwise),
    paste(sep = "\n",
      sprintf("stats <- dplyr::mutate(stats, .adj.p = p.adjust(.p.val, %s), .after = .p.val)", double_quote(p.adj)),
      "stats <- dplyr::arrange(stats, .p.val)" ))


  #________________________________________________________
  # Warnings raised while the tests run (including those the per-stratum
  # error handlers raise) are collected and reported once.
  #________________________________________________________
  msgs  <- character()
  stats <- withCallingHandlers(
    expr    = stats_run(df, regr, resp, stat.by, split.by, strata, pre, test, pairwise, p.adj),
    warning = function (w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart('muffleWarning') })

  stats_report(stats, unique(msgs), key_cols = c(split.by, stat.by))

  stats <- as_rbiom_tbl(stats)
  attr(stats, 'tbl_sum') <- c('Model' = model)
  attr(stats, 'code')    <- code

  # The data's provenance (its history and the command that made it), so
  # that `$code` is complete even when called directly; stats_table()
  # adds its own `cmd` afterwards.
  return (prov_attach(stats, data = df))
}



#________________________________________________________
# One summary of everything that went wrong while the tests ran:
# comparisons whose p-value is NA (their `key_cols` values are listed) and
# the unique messages raised by the underlying tests.
#________________________________________________________
stats_report <- function (stats, warnings, key_cols = NULL) {

  failed <- which(is.na(stats[['.p.val']]))

  if (length(failed) == 0 && length(warnings) == 0)
    return (invisible(NULL))

  labels <- character()
  if (length(failed) > 0) {
    key_cols <- intersect(key_cols, names(stats))
    labels   <- if (length(key_cols) == 0) "(all)" else
      do.call(paste, c(lapply(stats[failed, key_cols, drop = FALSE], as.character), sep = " / "))
  }

  # The messages are shown verbatim: escape braces so that cli doesn't
  # interpolate (and evaluate) anything in them, e.g. from group names.
  shown <- gsub("}", "}}", gsub("{", "{{", warnings, fixed = TRUE), fixed = TRUE)

  cli_warn(c(
    if (length(failed) > 0) c('!' = "Statistics could not be computed for {length(failed)} of {nrow(stats)} comparison{?s}: {.val {labels}}."),
    if (length(warnings) > 0) c('i' = "{qty(length(warnings))}Message{?s} from the underlying tests:"),
    stats::setNames(shown, rep('*', length(shown))) ))

  return (invisible(NULL))
}



#________________________________________________________
# Model formulas.
#
# Formulas are built as language objects (never parsed from text). The
# column names are the standardized ones (.resp, .regr, .stat.by, .strata)
# unless the caller gives others, e.g. corrplot_build()'s `y ~ x`.
#________________________________________________________

# Which fitting function the (fit, strata) combination uses.
stats_model_type <- function (fit, strata = FALSE) {
  if (fit == 'gam' || isTRUE(strata)) 'gam' else 'lm'
}


# The model formula. `regr`, `stat.by`, and `strata` are column names, or
# NULL to omit the term. `k` is the gam basis dimension (NULL to omit).
stats_formula <- function (fit, resp = ".resp", regr = ".regr", stat.by = NULL, strata = NULL, k = NULL) {

  stopifnot(fit %in% c('lm', 'log', 'gam'))

  smooth <- function (x, by = NULL) {
    args <- c(list(as.name('s'), as.name(x)), if (!is.null(by)) list(by = as.name(by)), list(bs = 'cs'), if (!is.null(k)) list(k = k))
    as.call(args)
  }

  rhs <- if (is.null(regr) && is.null(stat.by)) {
    1

  } else if (is.null(regr)) {
    as.name(stat.by)

  } else if (is.null(stat.by)) {
    switch(
      EXPR  = fit,
      'lm'  = as.name(regr),
      'log' = call('log', as.name(regr)),
      'gam' = smooth(regr) )

  } else {
    switch(
      EXPR  = fit,
      'lm'  = call('*', as.name(regr), as.name(stat.by)),
      'log' = call('*', call('log', as.name(regr)), as.name(stat.by)),
      'gam' = call('+', smooth(regr, by = stat.by), as.name(stat.by)) )
  }

  if (!is.null(strata))
    rhs <- call('+', rhs, as.call(list(as.name('s'), as.name(strata), bs = 're')))

  structure(call('~', as.name(resp), rhs), class = 'formula', .Environment = baseenv())
}


# The formula as code text, e.g. for the tibble header. `k` may be text
# such as "gam_k".
stats_formula_str <- function (fit, stat.by = FALSE, regr = TRUE, strata = FALSE, k = NULL) {

  frm <- stats_formula(
    fit     = fit,
    regr    = if (isTRUE(regr))    ".regr",
    stat.by = if (isTRUE(stat.by)) ".stat.by",
    strata  = if (isTRUE(strata))  ".strata",
    k       = if (!is.null(k)) as.name(k) )

  paste(deparse(frm, width.cutoff = 500L), collapse = "")
}


# Code that fits `model` from the data.frame named `var`.
stats_model_code <- function (var, fit, stat.by = FALSE, regr = TRUE, strata = FALSE) {

  type <- stats_model_type(fit, strata)

  if (type == 'lm')
    return (sprintf("model <- stats::lm(%s, data = %s)", stats_formula_str(fit, stat.by, regr, strata), var))

  paste(collapse = "\n", c(
    if (fit == 'gam') gam_basis_k_code(var, regr, stat.by),
    sprintf("model <- mgcv::gam(%s, data = %s, method = 'REML')", stats_formula_str(fit, stat.by, regr, strata, k = if (fit == 'gam') 'gam_k'), var) ))
}


#________________________________________________________
# Basis dimension for gam smooths: about one knot per five observations,
# at most 10, and fewer than the number of distinct x values in the
# smallest group (otherwise mgcv can't fit the smooth). Never below 3, the
# minimum for a cubic regression spline.
#________________________________________________________
gam_basis_k <- function (df, regr = ".regr", stat.by = NULL) {

  k <- min(10, nrow(df) / 5)

  if (!is.null(regr) && hasName(df, regr)) {
    x <- df[[regr]]
    u <- if (is.null(stat.by) || !hasName(df, stat.by)) { length(unique(x)) }
         else { min(lengths(lapply(split(x, df[[stat.by]], drop = TRUE), unique))) }
    k <- min(k, u - 1)
  }

  floor(max(3, k))
}

# The same computation as code text, for a data.frame named `var` with
# `.regr` / `.stat.by` columns.
gam_basis_k_code <- function (var, regr = TRUE, stat.by = FALSE) {

  if (!isTRUE(regr))
    return (sprintf("gam_k <- floor(max(3, min(10, nrow(%s) / 5)))", var))

  u <- if (isTRUE(stat.by)) {
    sprintf("min(lengths(lapply(split(%s$.regr, %s$.stat.by, drop = TRUE), unique)))", var, var)
  } else {
    sprintf("length(unique(%s$.regr))", var)
  }

  sprintf("gam_k <- floor(max(3, min(10, nrow(%s) / 5, %s - 1)))", var, u)
}


# The model description shown as the tibble header.
stats_model_str <- function (fit, stat.by = FALSE, regr = TRUE, strata = FALSE) {

  type <- stats_model_type(fit, strata)
  frm  <- stats_formula_str(fit, stat.by, regr, strata)

  if (type == 'lm') sprintf("lm(%s)", frm) else sprintf("gam(%s, method = 'REML')", frm)
}


# Goodness-of-fit columns, as code text and as values.
stats_gof_code <- function (type) {
  c(
    if (type == 'lm') c(.r.sqr = "gof$r.squared", .adj.r = "gof$adj.r.squared"),
    if (type == 'gam') c(.r.sqr = "gof$r.sq", .dev.expl = "gof$dev.expl"),
    .aic    = "stats::AIC(model)",
    .bic    = "stats::BIC(model)",
    .loglik = "as.numeric(stats::logLik(model))" )
}

stats_gof_vals <- function (type, model, gof) {
  c(
    if (type == 'lm') list(.r.sqr = gof$r.squared, .adj.r = gof$adj.r.squared),
    if (type == 'gam') list(.r.sqr = gof$r.sq, .dev.expl = gof$dev.expl),
    list(
      .aic    = stats::AIC(model),
      .bic    = stats::BIC(model),
      .loglik = as.numeric(stats::logLik(model)) ))
}



#________________________________________________________
# Fit the model used by stats_* and the corrplot trendlines. `strata` (a
# column name) adds a random intercept, fitted with mgcv::gam, and turns
# `lm`/`log` fits into linear mixed models. The data is embedded in the
# model's call so that emmeans can recover it.
#________________________________________________________
stats_fit_model <- function (df, fit, stat.by = FALSE, regr = ".regr", resp = ".resp", strata = NULL) {

  if (isTRUE(stat.by))  stat.by <- ".stat.by"
  if (isFALSE(stat.by)) stat.by <- NULL

  df <- as.data.frame(df)
  if (!is.null(strata)) df[[strata]] <- droplevels(as.factor(df[[strata]]))

  type <- stats_model_type(fit, !is.null(strata))
  k    <- if (fit == 'gam') gam_basis_k(df, regr, stat.by)
  frm  <- stats_formula(fit, resp, regr, stat.by, strata, k = k)

  if (type == 'lm') { do.call(stats::lm, list(frm, data = df)) }
  else              { do.call(mgcv::gam, list(frm, data = df, method = "REML")) }
}



#________________________________________________________
# Goodness-of-fit summary for an lm or gam model (base R; no broom).
#________________________________________________________
stats_glance <- function (model) {

  s <- summary(model)

  if (inherits(model, 'gam')) {

    tibble(
      adj.r.squared = s$r.sq,
      dev.expl      = s$dev.expl,
      edf           = sum(model$edf),
      logLik        = as.numeric(stats::logLik(model)),
      AIC           = stats::AIC(model),
      BIC           = stats::BIC(model),
      deviance      = stats::deviance(model),
      df.residual   = stats::df.residual(model),
      nobs          = stats::nobs(model) )

  } else {

    fs <- s$fstatistic

    tibble(
      r.squared     = s$r.squared,
      adj.r.squared = s$adj.r.squared,
      sigma         = s$sigma,
      statistic     = if (is.null(fs)) NA_real_ else unname(fs[['value']]),
      p.value       = if (is.null(fs)) NA_real_ else stats::pf(fs[['value']], fs[['numdf']], fs[['dendf']], lower.tail = FALSE),
      df            = if (is.null(fs)) NA_real_ else unname(fs[['numdf']]),
      logLik        = as.numeric(stats::logLik(model)),
      AIC           = stats::AIC(model),
      BIC           = stats::BIC(model),
      deviance      = stats::deviance(model),
      df.residual   = stats::df.residual(model),
      nobs          = stats::nobs(model) )
  }
}
