
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
    bdiv_boxplot   = list(biom = rare50, x = "Sex", stat.by = "Sex", permutations = 9),
    bdiv_corrplot  = list(biom = rare50, x = "Age", permutations = 9),
    bdiv_heatmap   = list(biom = rare50),
    bdiv_ord_plot  = list(biom = rare50, stat.by = "Sex", permutations = 9),
    bdiv_clusters  = list(biom = rare50),
    bdiv_distmat   = list(biom = rare50),
    bdiv_stats     = list(biom = rare50, stat.by = "Sex", permutations = 9),
    distmat_stats  = list(dm = bdiv_distmat(rare50), groups = pull(rare50, "Sex"), permutations = 9),
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

# Values a choice needs packages for.
schema_choice_skip <- function (fn, param, value) {
  (param == 'ord'    && value == 'tSNE' && !pkg_available('Rtsne')) ||
  (param == 'ord'    && value == 'UMAP' && !pkg_available('uwot'))  ||
  (param == 'format' && value == 'hdf5' && !pkg_available('h5lite')) ||
  (fn    == 'write_xlsx' && !pkg_available('openxlsx'))
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
  expect_true(length(schema) >= 23)

  for (fn in names(schema)) {

    entry <- schema[[fn]]
    fmls  <- formals(get(fn, envir = asNamespace('rbiom')))
    specs <- param_specs(entry)

    expect_identical(entry$fn, fn)
    expect_identical(names(specs), names(fmls), label = fn)
    expect_true(entry$returns$type %in% c("plot", "table", "factor", "matrix", "dist", "biom", "file"))

    for (nm in names(specs)) {

      spec  <- specs[[nm]]
      label <- paste0(fn, "(", nm, ")")

      expect_true(spec$type %in% types, label = label)
      expect_true(nzchar(spec$desc),    label = label)

      if (nm == '...') next

      if (identical(fmls[[nm]], quote(expr = ))) {
        expect_true(isTRUE(spec$required), label = label)
        expect_false(hasName(spec, 'default'), label = label)
        next
      }

      expect_true(hasName(spec, 'default'), label = label)

      if (is.null(spec$default_expr)) {
        expect_identical(spec$default, eval(fmls[[nm]], baseenv()), label = label)
      } else {
        expect_null(spec$default, label = label)
        expect_identical(spec$default_expr, paste(deparse(fmls[[nm]]), collapse = " "), label = label)
        expect_true(nzchar(spec$default_note %||% ""), label = label)
      }

      if (is.null(spec$default) && is.null(spec$default_expr))
        expect_true(hasName(spec, 'null'), label = label)

      for (cond in spec$when)
        expect_true(cond$param %in% names(specs), label = label)
    }
  }
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

  # Nothing but plain data.
  walk <- function (x) {
    expect_false(is.function(x) || is.environment(x) || is.language(x))
    if (is.list(x)) for (i in x) walk(i)
  }
  walk(schema)
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

      # On counts with zeros, Aitchison distances are NA without a
      # pseudocount (which only bdiv_distmat() takes), and Bhattacharyya
      # distances between samples with no taxa in common are Inf, which
      # clustering and ordination can't use.
      run_ids <- switch(
        EXPR          = fn,
        bdiv_distmat  = ids,
        bdiv_clusters = setdiff(ids, 'aitchison'),
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
      if (spec$type != 'field' || isTRUE(spec$unused)) next

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
  # `within` / `between` are not used.
  expect_error(fn_check("bdiv_stats", list(), d), "Either `stat.by`")
  expect_error(fn_check("bdiv_stats", list(stat.by = "Sex", test = "mantel"), d), "Choices are: adonis2, mrpp, betadisper")
  expect_identical(fn_check("bdiv_stats", list(regr = "==Age", test = "man"), d), list(regr = "Age", test = "mantel"))
  expect_identical(fn_check("bdiv_stats", list(stat.by = "!=Sex", within = "zzz"), d), list(stat.by = "Sex", within = "zzz"))

  # distmat_stats(): only the test options can be checked without R objects.
  expect_identical(fn_check("distmat_stats", list(test = "MRPP"), d), list(test = "mrpp"))
  expect_error(fn_check("distmat_stats", list(groups = 1:3), d), "R caller")
  expect_error(fn_check("distmat_stats", list(test = "zz"), d), "Invalid `test` option")

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

  # Every palette the schema suggests exists.
  for (entry in fn_schema())
    for (spec in entry$params)
      expect_true(all(spec$palettes %in% p$kind))
})
