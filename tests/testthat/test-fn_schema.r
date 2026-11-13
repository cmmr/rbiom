
# Arguments that make each covered function run quickly on the test data,
# and that every listed choice works with (e.g. `x` = `stat.by` so that
# every box plot test applies).
schema_base_args <- function (fn) {
  switch(
    EXPR           = fn,
    adiv_boxplot   = list(biom = rare50, x = "Sex", stat.by = "Sex"),
    adiv_corrplot  = list(biom = rare50, x = "Age"),
    adiv_matrix    = list(biom = rare50),
    adiv_stats     = list(biom = rare50, stat.by = "Sex", regr = "BMI"),
    adiv_table     = list(biom = rare50),
    bdiv_boxplot   = list(biom = rare50, x = "Sex", stat.by = "Sex", permutations = 9),
    bdiv_corrplot  = list(biom = rare50, x = "Age", permutations = 9),
    bdiv_heatmap   = list(biom = rare50),
    bdiv_ord_plot  = list(biom = rare50, stat.by = "Sex", permutations = 9),
    bdiv_ord_table = list(biom = rare50, stat.by = "Sex", permutations = 9),
    bdiv_clusters  = list(biom = rare50),
    bdiv_distmat   = list(biom = rare50),
    bdiv_stats     = list(biom = rare50, stat.by = "Sex", permutations = 9),
    bdiv_table     = list(biom = rare5),
    distmat_stats  = list(dm = bdiv_distmat(rare50), groups = pull(rare50, "Sex"), permutations = 9),
    distmat_ord_table = list(dm = bdiv_distmat(rare50)),
    rare_corrplot  = list(biom = hmp50),
    rare_multiplot = list(biom = hmp50),
    rare_stacked   = list(biom = hmp50),
    rarefy         = list(biom = hmp50),
    taxa_boxplot   = list(biom = rare50, x = "Sex", stat.by = "Sex", taxa = 3),
    taxa_corrplot  = list(biom = rare50, x = "Age", taxa = 2),
    taxa_heatmap   = list(biom = rare50),
    taxa_stacked   = list(biom = rare50),
    taxa_forest    = list(biom = hmp50, stat.by = "Sex"),
    taxa_da        = list(biom = hmp50, stat.by = "Sex"),
    taxa_clusters  = list(biom = rare50),
    taxa_matrix    = list(biom = rare50),
    taxa_stats     = list(biom = rare50, stat.by = "Sex", regr = "BMI", taxa = 2),
    taxa_table     = list(biom = rare50),
    write_biom     = list(biom = hmp5, file = tempfile(), overwrite = TRUE),
    write_xlsx     = list(biom = hmp5, file = tempfile(fileext = ".xlsx"), overwrite = TRUE),
    stop("No base arguments for ", fn) )
}

# Choices that other arguments must agree with.
schema_choice_args <- function (fn, param, value, args) {
  if (fn == 'bdiv_stats' && param == 'test' && value == 'mantel') {
    args$stat.by <- NULL
    args$regr    <- "Age"
  }
  if (fn == 'distmat_stats' && param == 'test' && value == 'mantel')
    args$groups <- pull(rare50, "Age")
  args
}

# Values that need a package that isn't installed, according to the
# schema's `needs` (the function's) and `choices_needs` (the value's).
schema_choice_skip <- function (fn, param, value) {
  entry <- fn_schema(fn)
  needs <- param_specs(entry)[[param]]$choices_needs
  pkgs  <- c(entry$needs, unlist(needs[intersect(value, names(needs))]))
  !all(vapply(pkgs, pkg_available, logical(1)))
}

schema_run <- function (fn, args) {
  suppressWarnings(suppressMessages(do.call(fn, args)))
}

# Like expect_no_error(), with a label saying which call failed.
expect_runs <- function (expr, label) {
  err <- tryCatch({ force(expr); NULL }, error = function (e) conditionMessage(e))
  expect(is.null(err), sprintf("%s failed: %s", label, err))
}

# The arguments fn_check() takes: all but the R-only ones.
check_args <- function (fn, args) {
  specs  <- param_specs(fn_schema(fn))
  r_only <- names(specs)[vapply(specs, function (s) isTRUE(s$r_only), logical(1))]
  args[setdiff(names(args), r_only)]
}

param_specs <- function (entry) {
  setNames(entry$params, vapply(entry$params, `[[`, character(1), 'name'))
}




test_that("each entry matches the function's formals", {

  types <- c(
    'biom', 'tree', 'dist', 'vector', 'file', 'dots', 'logical', 'integer', 'number', 'string',
    'choice', 'layers', 'field', 'rank', 'taxa', 'metric', 'aesthetic',
    'palette', 'tracks' )

  schema <- fn_schema()
  expect_true(length(schema) >= 32)

  # Every mismatch, as "fn(arg): problem". One expectation for the whole
  # schema: thousands of separate ones take half a minute.
  problems <- character(0)
  check    <- function (ok, label, problem) {
    if (!isTRUE(ok)) problems <<- c(problems, paste0(label, ": ", problem))
  }

  for (fn in names(schema)) {

    entry <- schema[[fn]]
    fmls  <- formals(get(fn, envir = asNamespace('rbiom')))
    specs <- param_specs(entry)

    check(identical(entry$fn, fn),               fn, "`fn` is wrong")
    check(identical(names(specs), names(fmls)),  fn, "arguments differ from formals()")
    check(isTRUE(entry$returns$type %in% c("plot", "table", "factor", "matrix", "dist", "biom", "file")), fn, "unknown return type")

    for (nm in names(specs)) {

      spec  <- specs[[nm]]
      label <- paste0(fn, "(", nm, ")")

      check(isTRUE(spec$type %in% types), label, "unknown type")
      check(isTRUE(nzchar(spec$desc)),    label, "no description")
      # A plain sentence: atima shows it as a tooltip and gives it to AI
      # models; the values are in `choices`.
      check(!grepl("`c(", spec$desc, fixed = TRUE), label, "description is R code, not a sentence")

      if (nm == '...') next

      if (identical(fmls[[nm]], quote(expr = ))) {
        check(isTRUE(spec$required),     label, "not marked as required")
        check(!hasName(spec, 'default'), label, "a required argument has a default")
        next
      }

      check(hasName(spec, 'default'), label, "no default")

      if (is.null(spec$default_expr)) {
        check(identical(spec$default, eval(fmls[[nm]], baseenv())), label, "default differs from formals()")
      } else {
        check(is.null(spec$default), label, "both `default` and `default_expr`")
        check(identical(spec$default_expr, paste(deparse(fmls[[nm]]), collapse = " ")), label, "`default_expr` differs from formals()")
        check(nzchar(spec$default_note %||% ""), label, "no `default_note`")
      }

      if (is.null(spec$default) && is.null(spec$default_expr))
        check(hasName(spec, 'null'), label, "a NULL default, but no `null`")

      for (cond in spec$when)
        check(isTRUE(cond$param %in% names(specs)), label, "`when` names an unknown argument")
    }
  }

  expect_identical(problems, character(0))
})




test_that("the schema survives a round trip through JSON", {

  inf_to_string <- function (x) {
    if (is.list(x)) return (lapply(x, inf_to_string))
    if (is.numeric(x) && any(is.infinite(x))) return (ifelse(x > 0, "Inf", "-Inf"))
    x
  }

  schema <- fn_schema()
  json   <- fn_schema(json = TRUE)

  rt <- jsonlite::fromJSON(json, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  expect_equal(rt, inf_to_string(schema))

  # List fields are arrays even with one value.
  raw <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  da  <- param_specs(raw$taxa_da)
  expect_type(da$method$choices, 'list')
  expect_type(param_specs(raw$adiv_boxplot)$colors$palettes, 'list')
  expect_type(param_specs(raw$bdiv_heatmap)$alpha$when[[1]]$has_any, 'list')

  # One function at a time.
  expect_identical(
    as.character(fn_schema("rarefy", json = TRUE)),
    as.character(jsonlite::toJSON(raw$rarefy, auto_unbox = TRUE, null = "null", digits = NA)) )

  # Nothing but plain data. (One expectation, not one per node.)
  not_data <- function (x) {
    if (is.function(x) || is.environment(x) || is.language(x)) return (TRUE)
    is.list(x) && any(vapply(x, not_data, logical(1)))
  }
  expect_false(not_data(schema))
})




test_that("every listed choice is accepted, and others are refused", {

  skip_on_cran()

  for (fn in names(fn_schema())) {

    specs <- param_specs(fn_schema(fn))
    base  <- schema_base_args(fn)

    for (nm in names(specs)) {

      spec <- specs[[nm]]
      if (isTRUE(spec$r_only) || is.null(spec$choices)) next
      if (spec$type == 'metric') next   # see below

      values <- if (spec$type == 'layers') {
        c(vapply(spec$choices, `[[`, character(1), 'name'),
          paste(vapply(spec$choices, `[[`, character(1), 'letter'), collapse = ""))
      } else {
        unlist(spec$choices)
      }

      for (value in values) {
        if (schema_choice_skip(fn, nm, value)) next
        args <- schema_choice_args(fn, nm, value, base)
        args[[nm]] <- value
        label <- sprintf("%s(%s = %s)", fn, nm, value)
        expect_runs(schema_run(fn, args), label)
        expect_runs(fn_check(fn, check_args(fn, args), args$biom %||% rare50), paste("fn_check:", label))
      }

      if (schema_choice_skip(fn, nm, "")) next
      args <- base
      args[[nm]] <- "zzzz"
      expect_error(schema_run(fn, args), label = sprintf("%s(%s = 'zzzz')", fn, nm))
      expect_error(fn_check(fn, check_args(fn, args), args$biom %||% rare50))
    }
  }
})




test_that("every metric is accepted, and others are refused", {

  skip_on_cran()

  for (fn in names(fn_schema())) {

    specs <- param_specs(fn_schema(fn))
    base  <- schema_base_args(fn)

    for (nm in names(specs)) {

      spec <- specs[[nm]]
      if (spec$type != 'metric') next

      ids  <- vapply(spec$choices, `[[`, character(1), 'id')
      args <- base[intersect(names(base), names(specs))]

      # Skip the statistics where that's an option.
      if ('none' %in% unlist(specs$test$choices)) args$test <- 'none'

      # Aitchison distances are NA (ecodive 2.3.1 gives NaN for every
      # pair), and Bhattacharyya distances between samples with no taxa
      # in common are Inf. The figures, ordinations, clusterings, and
      # tests refuse both (see "distances that can't be used").
      run_ids <- switch(
        EXPR          = fn,
        bdiv_distmat  = ids,
        setdiff(ids, c('aitchison', 'bhattacharyya')) )

      # One call per metric where only one is allowed; otherwise one call
      # with all of them.
      if (identical(spec$max_n, 1)) {
        for (id in run_ids) {
          args[[nm]] <- id
          expect_runs(schema_run(fn, args), sprintf("%s(%s = %s)", fn, nm, id))
        }
      } else {
        args[[nm]] <- run_ids
        expect_runs(schema_run(fn, args), sprintf("%s(%s = <all>)", fn, nm))
      }

      d <- fn_check_data(rare50)
      for (id in ids) {
        args[[nm]] <- id
        expect_runs(fn_check(fn, check_args(fn, args), d), sprintf("fn_check: %s(%s = %s)", fn, nm, id))
      }

      args[[nm]] <- "zzzz"
      expect_error(schema_run(fn, args))
      expect_error(fn_check(fn, check_args(fn, args), d))
    }
  }
})




test_that("fields: fn_check() follows `class` and `coerce`, as the functions do", {

  d    <- jsonlite::fromJSON(jsonlite::toJSON(fn_check_data(hmp50), auto_unbox = TRUE))
  cats <- "Body Site"
  nums <- "Age"

  for (fn in names(fn_schema())) {

    specs <- param_specs(fn_schema(fn))
    base  <- schema_base_args(fn)

    for (nm in names(specs)) {

      spec <- specs[[nm]]
      if (spec$type != 'field' || !is.null(spec$deprecated)) next

      # Other arguments that the field must not clash with.
      args <- base[setdiff(names(base), c('biom', 'file', 'overwrite', nm))]
      if (nm %in% c('x', 'stat.by', 'facet.by', 'strata', 'within', 'between') && fn %in% c('adiv_boxplot', 'bdiv_boxplot', 'taxa_boxplot'))
        args$x <- args$stat.by <- NULL
      if (nm == 'regr') args$stat.by <- NULL

      for (field in c(cats, nums)) {

        ok <- spec$class == 'any' || isTRUE(spec$coerce) ||
          (spec$class == 'categorical' && field %in% cats) ||
          (spec$class == 'numeric'     && field %in% nums)

        a <- args
        a[[nm]] <- field
        label <- sprintf("fn_check('%s', list(%s = '%s'))", fn, nm, field)

        if (ok) { expect_runs(fn_check(fn, a, d), label) }
        else    { expect_error(fn_check(fn, a, d), label = label) }
      }

      a <- args
      a[[nm]] <- "zzzz"
      expect_error(fn_check(fn, a, d), "zzzz")
    }
  }


  # Numeric fields where categorical ones are expected, when `coerce`.
  skip_on_cran()
  expect_no_error(schema_run('adiv_boxplot',  list(rare50, x = "Age", stat.by = "Age", facet.by = "BMI")))
  expect_no_error(schema_run('taxa_boxplot',  list(rare50, x = "Age", taxa = 2)))
  expect_no_error(schema_run('bdiv_boxplot',  list(rare50, x = "Age", test = "none")))
  expect_no_error(schema_run('adiv_corrplot', list(rare50, x = "Age", stat.by = "BMI", strata = "BMI")))
  expect_no_error(schema_run('bdiv_corrplot', list(rare50, x = "Age", stat.by = "BMI", permutations = 9)))
  expect_no_error(schema_run('taxa_corrplot', list(rare50, x = "Age", stat.by = "BMI", taxa = 2)))
  expect_no_error(schema_run('taxa_stacked',  list(rare50, label.by = "Age", order.by = "BMI")))
  expect_no_error(schema_run('taxa_da',       list(hmp50, stat.by = "Sex", covariates = "Age")))
  expect_no_error(schema_run('taxa_heatmap',  list(rare50, tracks = "Age")))
  expect_no_error(schema_run('taxa_stats',    list(rare50, stat.by = ".taxa", split.by = "Sex", test = "kruskal", taxa = 3)))
})




test_that("`when` conditions: an argument has no effect while they fail", {

  skip_on_cran()

  # A value that fails the condition.
  falsify <- function (cond, specs) {
    spec <- specs[[cond$param]]
    if (!is.null(cond[['not_null']])) return (list(NULL))
    if (!is.null(cond[['is_null']]))  return (list("Sex"))
    if (!is.null(cond[['is_not']]))   return (list(cond[['is_not']][[1]]))
    choices <- if (spec$type == 'layers') {
      vapply(spec$choices, `[[`, character(1), 'name')
    } else if (spec$type == 'metric') {
      c("bray", "shannon")
    } else {
      unlist(spec$choices)
    }
    list(setdiff(choices, c(cond[['is']], cond[['has_any']]))[[1]])
  }

  # A value other than the default (and than the base arguments').
  other_value <- function (spec) {
    switch(
      EXPR = spec$type,
      aesthetic = !isTRUE(spec[['default']]),
      logical   = TRUE,
      integer   = 99,
      field     = if (identical(spec$class, 'numeric')) "Age" else "Body Site",
      number    = if (is.null(spec[['default']])) 30
                  else if (is.null(spec[['max']])) spec[['default']] * 2 + 1
                  else spec[['min']] + 0.9 * (spec[['max']] - spec[['min']]),
      choice    = setdiff(unlist(spec$choices), spec[['default']])[[1]],
      taxa      = 2 )
  }

  comparable <- function (x) {
    if (inherits(x, 'ggplot')) return (list(
      data   = as.data.frame(x$data),
      stats  = if (!is.null(x$stats)) as.data.frame(x$stats),
      layers = lapply(ggplot2::ggplot_build(x)$data, as.data.frame) ))
    if (inherits(x, 'dist')) return (as.matrix(x))
    if (is.matrix(x)) return (unclass(as.matrix(x)))
    x
  }

  n <- 0
  for (fn in names(fn_schema())) {

    specs <- param_specs(fn_schema(fn))

    for (nm in names(specs)) {

      spec <- specs[[nm]]
      if (is.null(spec$when)) next

      cond <- spec$when[[1]]
      args <- schema_base_args(fn)
      args[cond$param] <- falsify(cond, specs)

      a <- args
      b <- args
      b[nm] <- list(other_value(spec))

      label <- sprintf("%s(%s) with %s = %s", fn, nm, cond$param, deparse(args[[cond$param]]))
      expect_equal(comparable(schema_run(fn, a)), comparable(schema_run(fn, b)), ignore_attr = TRUE, label = label)
      n <- n + 1
    }
  }
  expect_true(n > 10)
})




test_that("fn_check()", {

  d <- fn_check_data(hmp50)

  expect_named(d, c('fields', 'ranks', 'tree'))
  expect_identical(d$fields$Age$class, "numeric")
  expect_identical(d$fields$Sex$class, "categorical")
  expect_true(d$tree)

  # Conversions.
  x <- fn_check("adiv_boxplot", list(x = "body", adiv = c("shan", "Simpson"), layers = "xd"), hmp50)
  expect_identical(x, list(x = "Body Site", adiv = c("shannon", "simpson"), layers = c("box", "dot")))

  x <- fn_check("taxa_boxplot", list(rank = 2, taxa = 3), d)
  expect_identical(x, list(rank = "Phylum", taxa = 3))

  x <- fn_check("bdiv_boxplot", list(x = "==Body Site", bdiv = "Bray"), d)
  expect_identical(x, list(x = "Body Site", bdiv = "bray", within = "Body Site"))

  x <- fn_check("taxa_heatmap", list(clust = "ward"), d)
  expect_identical(x, list(clust = "ward.D"))

  x <- fn_check("write_biom", list(format = "TAB"), d)
  expect_identical(x, list(format = "tab"))

  # The function's own messages.
  expect_error(fn_check("adiv_boxplot", list(ci = "zz"), d), "Invalid `ci` option")
  expect_error(fn_check("adiv_boxplot", list(x = "Age", strata = "nope"), d), "Invalid `strata` option")
  expect_error(fn_check("rare_corrplot", list(stat.by = "Age"), d), "not categorical")
  expect_error(fn_check("adiv_corrplot", list(x = "Sex"), d), "not numeric")
  expect_error(fn_check("adiv_corrplot", list(x = "Age", level = 2), d), "between 0.5 and 1")
  expect_error(fn_check("bdiv_boxplot", list(permutations = 1.5), d), "whole numbers")
  expect_error(fn_check("rare_stacked", list(rline = "a"), d), "TRUE, FALSE, or an integer")
  expect_error(fn_check("taxa_da", list(), d), "Give either `stat.by`")
  expect_error(fn_check("taxa_da", list(stat.by = "Sex", regr = "Age"), d), "not both")
  expect_error(fn_check("write_biom", list(format = "xml"), d), "should be one of")

  # Statistics: arguments a test needs, or can't use.
  expect_error(fn_check("adiv_stats", list(test = "kruskal"), d), "`stat.by` is required")
  expect_error(fn_check("adiv_stats", list(test = "emtrends"), d), "`regr` is required")
  expect_error(fn_check("adiv_stats", list(stat.by = "Sex", test = "wilcox", strata = "Body Site"), d), "not supported by rank tests")
  expect_error(fn_check("taxa_stats", list(stat.by = "Sex", test = "kruskal", covariates = "Age"), d), "not supported by rank tests")
  expect_error(adiv_stats(rare50, test = "emtrends"), "`regr` is required")

  # ".taxa" compares taxa in taxa_stats(), but isn't a field elsewhere.
  x <- fn_check("taxa_stats", list(stat.by = ".taxa", split.by = c("sex", ".taxa")), d)
  expect_identical(x, list(stat.by = ".taxa", split.by = c("Sex", ".taxa")))
  expect_error(fn_check("taxa_stats", list(stat.by = "Sex", strata = ".taxa"), d), "Invalid `strata` option")
  expect_error(fn_check("taxa_boxplot", list(x = ".taxa"), d), "Invalid `x` option")

  # bdiv_stats(): the tests depend on `regr`; prefixes are dropped, and
  # `within` / `between` are deprecated (refused).
  expect_error(fn_check("bdiv_stats", list(), d), "Either `stat.by`")
  expect_error(fn_check("bdiv_stats", list(stat.by = "Sex", test = "mantel"), d), "Choices are: adonis2, mrpp, betadisper")
  expect_identical(fn_check("bdiv_stats", list(regr = "==Age", test = "man"), d), list(regr = "Age", test = "mantel"))
  expect_identical(fn_check("bdiv_stats", list(stat.by = "!=Sex"), d), list(stat.by = "Sex"))
  expect_error(fn_check("bdiv_stats", list(stat.by = "Sex", within = "Sex"), d), "deprecated")

  # distmat_stats(): only the test options can be checked without R objects.
  expect_identical(fn_check("distmat_stats", list(test = "MRPP"), d), list(test = "mrpp"))
  expect_error(fn_check("distmat_stats", list(groups = 1:3), d), "R caller")
  expect_error(fn_check("distmat_stats", list(test = "zz"), d), "Invalid `test` option")

  # Tables: ".all" adds every field, but only on its own.
  expect_setequal(fn_check("adiv_table", list(md = ".all"), d)$md, setdiff(names(d$fields), ".sample"))
  expect_error(fn_check("adiv_table", list(md = c(".all", "Sex")), d), "Invalid `md` option")
  expect_error(adiv_table(rare5, md = c(".all", "Sex")), "Invalid `md` option")
  expect_identical(
    fn_check("bdiv_table", list(md = c("==Body Site", "Age"), delta = "age"), d),
    list(md = c("Body Site", "Age"), delta = "Age", within = "Body Site") )
  expect_identical(fn_check("bdiv_ord_table", list(ord = "nm", k = 3), d), list(ord = "NMDS", k = 3))
  expect_error(fn_check("bdiv_ord_table", list(k = 4), d), "between 2 and 3")
  expect_error(fn_check("distmat_ord_table", list(dm = 1), d), "R caller")

  # A tree is needed for phylogenetic metrics.
  d0 <- d
  d0$tree <- FALSE
  expect_error(fn_check("adiv_boxplot", list(adiv = "faith"), d0), "phylogenetic tree")
  expect_error(fn_check("bdiv_heatmap", list(bdiv = c("bray", "unweighted_unifrac")), d0), "phylogenetic tree")
  expect_no_error(fn_check("adiv_boxplot", list(adiv = "faith"), d))
  expect_error(adiv_boxplot(min50, adiv = "faith"), "phylogenetic tree")

  # The request itself.
  expect_error(fn_check("adiv_boxplot", list(zz = 1), d), "no argument")
  expect_error(fn_check("bdiv_boxplot", list(cpus = 2), d), "R caller")
  expect_error(fn_check("adiv_boxplot", list(1), d), "named")
  expect_error(fn_check("adiv_corrplot", list(), d), "needs")
  expect_error(fn_check("stats_table", list(), d), "not described")
  expect_error(fn_check("adiv_boxplot", list(), list()), "fn_check_data")
})




test_that("fixes that came with the schema", {

  skip_on_cran()

  # Ward clustering: "ward.D" and "ward.D2", and "ward" for "ward.D".
  expect_no_error(taxa_heatmap(rare50, clust = "ward.D2"))
  expect_no_message(p <- taxa_heatmap(rare50, clust = "ward"))
  expect_match(paste(p$code, collapse = "\n"), "ward.D", fixed = TRUE)
  expect_no_error(taxa_stacked(rare50, clust = "ward.D"))

  # Tests that don't apply to the plotted values are refused.
  expect_error(adiv_boxplot(rare50, x = "Sex", test = "adonis2"), "Invalid `test` option")
  expect_error(taxa_boxplot(rare50, x = "Sex", test = "mrpp"),    "Invalid `test` option")
  expect_error(bdiv_boxplot(rare50, x = "Sex", test = "wilcox"),  "Invalid `test` option")
  expect_error(adiv_corrplot(rare50, x = "Age", test = "mantel"), "Invalid `test` option")
  expect_error(rare_corrplot(hmp50, test = "mantel"),             "Invalid `test` option")

  # taxa_stacked()'s xlab.angle default is a string, like the choices.
  expect_identical(formals(taxa_stacked)$xlab.angle, "90")
  expect_no_error(taxa_stacked(rare50, xlab.angle = 30))

  # Heat map options give the usual message.
  expect_error(taxa_heatmap(rare50, legend = "top"), "Invalid `legend` option")
  expect_error(bdiv_heatmap(rare50, rescale = "z"),  "Invalid `rescale` option")
  expect_error(bdiv_heatmap(rare50, label = "z"),    "Invalid `label` option")
  expect_error(bdiv_heatmap(rare50, title = FALSE),  "must be TRUE, NULL")
})




test_that("fn_check() refuses deprecated arguments, saying what to use instead", {

  d <- fn_check_data(rare50)
  n <- 0

  for (fn in names(fn_schema())) {
    for (spec in fn_schema(fn)$params) {

      if (is.null(spec$deprecated)) next
      n <- n + 1

      args  <- setNames(list(if (spec$type == 'logical') TRUE else "Sex"), spec$name)
      label <- sprintf("fn_check(%s, %s)", fn, spec$name)
      err   <- tryCatch(fn_check(fn, args, d), error = function (e) conditionMessage(e))

      expect(is.character(err), paste(label, "was accepted"))
      if (!is.character(err)) next
      expect_match(err, paste0("`", spec$name, "` is deprecated"), fixed = TRUE, label = label)
      expect_match(err, substr(spec$note, 1, 30), fixed = TRUE, label = label)
    }
  }

  expect_gte(n, 12)

  err <- tryCatch(fn_check("bdiv_boxplot", list(weighted = TRUE), d), error = function (e) conditionMessage(e))
  expect_match(err, "since rbiom 2.3.0", fixed = TRUE)
  expect_match(err, "Name the metric in `bdiv`", fixed = TRUE)

  # The functions themselves still accept them, with a warning.
  expect_warning(bdiv_distmat(rare50, bdiv = "unifrac", weighted = TRUE), "deprecated")
})




test_that("a bad metric is reported under the caller's argument", {

  d <- fn_check_data(rare50)

  for (fn in c("adiv_boxplot", "adiv_stats", "adiv_table")) {
    expect_error(fn_check(fn, list(adiv = "nosuch"), d), '`adiv = "nosuch"` does not match any option', fixed = TRUE)
    expect_error(fn_check(fn, list(adiv = "s"),      d), '`adiv = "s"` matches multiple options', fixed = TRUE)
  }
  for (fn in c("bdiv_boxplot", "bdiv_distmat", "bdiv_ord_plot")) {
    expect_error(fn_check(fn, list(bdiv = "nosuch"), d), '`bdiv = "nosuch"` does not match any option', fixed = TRUE)
  }

  # The functions give the same message.
  expect_error(adiv_boxplot(rare50, adiv = "nosuch"), '`adiv = "nosuch"` does not match any option', fixed = TRUE)
  expect_error(adiv_matrix(rare50,  adiv = "nosuch"), '`adiv = "nosuch"` does not match any option', fixed = TRUE)
  expect_error(bdiv_boxplot(rare50, bdiv = "nosuch"), '`bdiv = "nosuch"` does not match any option', fixed = TRUE)
  expect_error(bdiv_distmat(rare50, bdiv = "nosuch"), '`bdiv = "nosuch"` does not match any option', fixed = TRUE)
  expect_error(bdiv_distmat(rare50, bdiv = "1"),      '`bdiv` must contain alphabetic characters', fixed = TRUE)

  err <- tryCatch(bdiv_heatmap(rare50, bdiv = "nosuch"), error = function (e) conditionMessage(e))
  expect_no_match(err, "metric =", fixed = TRUE)
  expect_match(err, "bray", fixed = TRUE)
})




test_that("choices and functions that need an optional package are marked", {

  schema <- fn_schema()

  expect_identical(param_specs(schema$bdiv_ord_plot)$ord$choices_needs,     list(tSNE = "Rtsne", UMAP = "uwot"))
  expect_identical(param_specs(schema$bdiv_ord_table)$ord$choices_needs,    list(tSNE = "Rtsne", UMAP = "uwot"))
  expect_identical(param_specs(schema$distmat_ord_table)$ord$choices_needs, list(tSNE = "Rtsne", UMAP = "uwot"))
  expect_identical(param_specs(schema$bdiv_ord_plot)$layers$choices_needs,  list(density = "ggdensity"))
  expect_identical(param_specs(schema$write_biom)$format$choices_needs,     list(hdf5 = "h5lite"))
  expect_identical(param_specs(schema$taxa_da)$strata$needs,                "lmerTest")
  expect_identical(param_specs(schema$taxa_forest)$strata$needs,            "lmerTest")
  expect_identical(schema$write_xlsx$needs,                                 "openxlsx")

  # `choices` stay plain strings (or letter/name pairs for layers); every
  # marked choice is one of them, and every package is a suggested one.
  desc     <- read.dcf(system.file("DESCRIPTION", package = "rbiom"), fields = "Suggests")
  suggests <- trimws(sub("\\(.*", "", strsplit(desc[[1]], ",")[[1]]))
  pkgs     <- character(0)

  for (entry in schema) {
    pkgs <- c(pkgs, entry$needs)
    for (spec in entry$params) {
      pkgs <- c(pkgs, spec$needs)
      if (is.null(cn <- spec$choices_needs)) next
      names <- if (spec$type == 'layers') vapply(spec$choices, `[[`, character(1), 'name') else unlist(spec$choices)
      expect(all(names(cn) %in% names), sprintf("%s(%s): `choices_needs` names unknown choices", entry$fn, spec$name))
      pkgs <- c(pkgs, unlist(cn))
    }
  }
  expect_true(all(pkgs %in% suggests))

  # An object in JSON.
  raw <- jsonlite::fromJSON(fn_schema("bdiv_ord_plot", json = TRUE), simplifyVector = FALSE)
  expect_identical(param_specs(raw)$ord$choices_needs, list(tSNE = "Rtsne", UMAP = "uwot"))
  expect_identical(param_specs(raw)$ord$choices, list("PCoA", "tSNE", "NMDS", "UMAP"))
})




test_that("entries say which counts the function wants", {

  counts <- vapply(fn_schema(), function (e) e$counts %||% NA_character_, character(1))

  expect_true(all(counts %in% c("rarefied", "unrarefied", NA)))

  expect_setequal(
    names(which(counts == "unrarefied")),
    c("rare_corrplot", "rare_multiplot", "rare_stacked", "rarefy", "taxa_forest", "taxa_da", "write_xlsx") )

  expect_setequal(names(which(is.na(counts))), c("write_biom", "distmat_stats", "distmat_ord_table"))

  for (fn in c("adiv_boxplot", "adiv_stats", "bdiv_ord_plot", "bdiv_distmat", "taxa_boxplot", "taxa_heatmap", "taxa_stats"))
    expect_identical(counts[[fn]], "rarefied", label = fn)

  # taxa_da() says so for rarefied counts, and not for unrarefied ones.
  expect_message(taxa_da(rare50, stat.by = "Sex"), "unrarefied counts are better")
  expect_no_message(taxa_da(hmp50, stat.by = "Sex"), message = "unrarefied counts are better")
})




test_that("distances that can't be used are refused with a clear message", {

  skip_on_cran()

  # Bhattacharyya: Inf between samples with no taxa in common.
  dm <- bdiv_distmat(rare50, bdiv = "bhattacharyya")
  expect_true(any(is.infinite(dm)))

  msg <- function (expr) tryCatch({ force(expr); "" }, error = function (e) conditionMessage(e))

  x <- msg(bdiv_heatmap(rare50, bdiv = "bhattacharyya"))
  expect_match(x, '"bhattacharyya" metric gave infinite distances for 279 of 1176 sample pairs', fixed = TRUE)
  expect_match(x, "can't be used for a heat map", fixed = TRUE)
  expect_match(x, "no taxa in common", fixed = TRUE)

  expect_match(msg(bdiv_ord_plot(rare50, bdiv = "bhattacharyya")),                 "can't be used for an ordination",  fixed = TRUE)
  expect_match(msg(bdiv_ord_table(rare50, bdiv = "bhattacharyya")),                "can't be used for an ordination",  fixed = TRUE)
  expect_match(msg(bdiv_clusters(rare50, bdiv = "bhattacharyya")),                 "can't be used for clustering",     fixed = TRUE)
  expect_match(msg(bdiv_boxplot(rare50, bdiv = "bhattacharyya", x = "Sex")),       "can't be used for a box plot",     fixed = TRUE)
  expect_match(msg(bdiv_corrplot(rare50, bdiv = "bhattacharyya", x = "Age")),      "can't be used for a scatter plot", fixed = TRUE)
  expect_match(msg(bdiv_stats(rare50, bdiv = "bhattacharyya", stat.by = "Sex")),   "can't be used for statistics",     fixed = TRUE)

  # Several metrics: the one that can't be used is named.
  expect_match(msg(bdiv_boxplot(rare50, bdiv = c("bray", "bhattacharyya"), x = "Sex")), '"bhattacharyya" metric', fixed = TRUE)

  # Missing distances are named as such.
  expect_match(
    object = tryCatch(validate_finite_dist(c(1, NA, 2), "x", "a test"), error = function (e) conditionMessage(e)),
    regexp = '"x" metric gave missing (NA) distances for 1 of 3 sample pairs', fixed = TRUE )
  expect_no_error(validate_finite_dist(c(0.1, -Inf, 2), "x", "a test"))

  # On samples that share taxa, the metric works.
  b <- rarefy(hmp50[pull(hmp50, "Body Site") == "Stool"])
  if (all(is.finite(bdiv_distmat(b, bdiv = "bhattacharyya"))))
    expect_no_error(bdiv_heatmap(b, bdiv = "bhattacharyya"))
})




test_that("missing distances are found where bdiv_table() drops them", {

  skip_on_cran()

  # A metric that gives one NA distance. (No cache: the mock must run.)
  op <- options(rbiom.cache_dir = FALSE)
  on.exit(options(op))
  local_mocked_bindings(bdiv_distmat = function (biom, bdiv, ...) {
    dm <- ecodive::bray(as.matrix(biom$counts), margin = 2L)
    dm[2] <- NA
    dm })

  msg <- function (expr) tryCatch({ force(expr); "" }, error = function (e) conditionMessage(e))

  expect_match(msg(bdiv_boxplot(rare50,  bdiv = "bray", x = "Sex")), '"bray" metric gave missing (NA) distances for 1 of 1176', fixed = TRUE)
  expect_match(msg(bdiv_corrplot(rare50, bdiv = "bray", x = "Age")), '"bray" metric gave missing (NA) distances for 1 of 1176', fixed = TRUE)

  # With `within` / `between`, the NA rows are the pairs left out.
  expect_no_error(suppressWarnings(bdiv_boxplot(rare50, bdiv = "bray", x = "==Sex", test = "none")))
})




test_that("list_palettes()", {

  p <- list_palettes()
  expect_s3_class(p, 'tbl_df')
  expect_named(p, c('name', 'kind', 'colors'))
  expect_setequal(p$name, names(PALETTES))
  expect_false(anyNA(p$kind))
  expect_true(all(grepl("^#[0-9A-F]{6}$", unlist(p$colors))))

  expect_setequal(list_palettes("diverging")$kind, "diverging")
  expect_setequal(list_palettes(c("seq", "disc"))$kind, c("sequential", "discrete"))
  expect_error(list_palettes("zz"))

  # Every palette kind the schema suggests exists.
  kinds <- unlist(lapply(fn_schema(), function (e) lapply(e$params, `[[`, 'palettes')))
  expect_true(length(kinds) > 0)
  expect_true(all(kinds %in% p$kind))
})
