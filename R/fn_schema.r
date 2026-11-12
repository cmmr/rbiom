#' Describe the arguments of rbiom functions.
#'
#' A machine-readable description of what each argument of an rbiom
#' function accepts: its kind, its fixed choices, numeric limits, whether
#' `NULL` is allowed and what it means, and where dataset-dependent
#' choices come from. Another program (such as a web interface) can use it
#' to draw input controls and to check a request without calling the
#' function. The same description is what the functions use to check their
#' own arguments, so the two cannot disagree.
#'
#' @param fn   The name of a function, e.g. `"adiv_boxplot"`. `NULL`
#'        describes every covered function. Default: `NULL`
#'
#' @param json   Return the description as JSON text instead of an R list.
#'        Default: `FALSE`
#'
#' @return With `fn`, one entry (a list). Without, a named list of entries,
#'         one per covered function, so `names(fn_schema())` lists the
#'         covered functions. With `json = TRUE`, the same as a JSON string.
#'
#'
#' @section Entries:
#'
#' Each entry has:
#' \describe{
#'   \item{`fn`}{ The function's name. }
#'   \item{`returns`}{ What the function returns: `type` is one of `"plot"`,
#'     `"table"`, `"factor"`, `"matrix"`, `"dist"`, `"biom"`, or `"file"`;
#'     `has` lists the `$` fields the result carries (e.g. `"data"`,
#'     `"code"`, `"stats"`). }
#'   \item{`rules`}{ (Optional) Constraints between arguments, each one of:
#'     `one_of` (exactly one of these arguments must be non-`NULL`);
#'     `any_of` (at least one must be); or `when` (conditions, as for an
#'     argument's `when`) with `requires` (these arguments must then be
#'     non-`NULL`) or `forbids` (they must then be `NULL`). }
#'   \item{`params`}{ One entry per argument, in the order of [formals()]. }
#' }
#'
#'
#' @section Arguments:
#'
#' Every argument has `name`, `type`, and `desc` (the first sentence of
#' its help text). The other fields are present only when they apply.
#'
#' \describe{
#'   \item{`type`}{ One of:
#'     `"logical"`, `"integer"`, `"number"`, `"string"`;
#'     `"choice"` (one or more of `choices`);
#'     `"layers"` (plot layers: `choices` lists each `letter` and `name`;
#'       letters may be run together, as in `"xd"`);
#'     `"field"` (metadata field names: see `class`);
#'     `"rank"` (taxonomic ranks of the dataset, by name or by position:
#'       `1` is the highest rank, `-1` the lowest, `0` the OTUs);
#'     `"taxa"` (see `forms`);
#'     `"metric"` (`div` is `"alpha"` or `"beta"`; `choices` lists each
#'       metric's `id`, `name`, whether it is `phylo`genetic, needing a
#'       tree, and whether it needs integer counts, `int_only`);
#'     `"aesthetic"` (colors, shapes, or patterns: `aes`, and `by`, the
#'       argument whose values they are mapped to);
#'     `"palette"`; `"tracks"`;
#'     `"biom"`, `"tree"`, `"dist"`, `"vector"`, `"file"`, `"dots"` (R
#'       objects, a file path, and `...`: see `r_only`). }
#'   \item{`required`}{ `TRUE` when the argument has no default. }
#'   \item{`default`}{ The default value, from [formals()]. }
#'   \item{`default_expr`, `default_note`}{ For a default that is an
#'     expression, such as `stat.by = x`: the expression (not evaluated)
#'     and what it means. `default` is then `NULL`. }
#'   \item{`null`}{ Present when `NULL` is accepted: what it means. }
#'   \item{`choices`}{ The fixed options. Matching ignores case and
#'     accepts unambiguous abbreviations. }
#'   \item{`choices_when`}{ Options that depend on other arguments: a list
#'     of `when` conditions (as below) with the `choices` that apply while
#'     they hold. The first that holds applies; otherwise `choices`. }
#'   \item{`min`, `max`}{ Numeric limits. `"integer"` values must be whole. }
#'   \item{`max_n`}{ How many values may be given. Absent: any number. }
#'   \item{`per_side`}{ One value for both rows and columns of a heat
#'     map, or two: for rows, then columns. }
#'   \item{`class`}{ For fields: `"categorical"`, `"numeric"`, or `"any"`. }
#'   \item{`coerce`}{ For categorical fields: numeric fields are accepted
#'     too, and treated as categories. }
#'   \item{`ops`}{ For fields: field names may be prefixed with `"=="`
#'     (compare within groups) or `"!="` (between groups). }
#'   \item{`extra`}{ For fields: values accepted besides the dataset's
#'     fields, such as `".taxa"`. }
#'   \item{`forms`}{ The forms accepted, when there are several:
#'     `"logical"` (`TRUE` or `FALSE`); `"true"` (`TRUE` only);
#'     `"integer"`; `"string"`; `"choice"` (one of `choices`);
#'     `"palette"` (a palette name from [list_palettes()], optionally
#'       prefixed with `"-"` to reverse it);
#'     `"values"` (a vector of colors, shapes, or patterns);
#'     `"named_values"` (a vector named by the values of `by`);
#'     `"count"` (the n most abundant taxa); `"fraction"` (a number
#'       below 1: taxa with at least that mean relative abundance);
#'     `"names"` (taxa names at the chosen rank);
#'     `"fields"` (metadata field names); `"spec"` (a list of options;
#'       see the function's help). }
#'   \item{`palettes`}{ The palette kinds suited to the argument. }
#'   \item{`needs`}{ An optional package the argument requires. }
#'   \item{`when`}{ Conditions under which the argument has any effect;
#'     all must hold. Each names a `param` and one of: `is` (its value is
#'     one of these), `is_not`, `has_any` (it includes one of these),
#'     `is_null`, `not_null`. Values are canonical: layer names, not
#'     letters, and metric ids. Absent: the argument may always matter. }
#'   \item{`deprecated`}{ The version that deprecated the argument. }
#'   \item{`unused`}{ Accepted for compatibility, but has no effect. }
#'   \item{`r_only`}{ Only an R caller can set it: an R object, a file
#'     path, a host setting such as `cpus`, or `...`. [fn_check()] refuses
#'     it. }
#'   \item{`note`}{ Anything else a caller needs to know. }
#' }
#'
#' In JSON, `NULL` is `null`, numbers are written in full, `Inf` is the
#' string `"Inf"`, and the list fields (`choices`, `forms`, `ops`,
#' `palettes`, `has`, and the values in `when`) are always arrays.
#'
#'
#' @seealso [fn_check()] to check arguments against a dataset;
#'   [list_palettes()] and [list_metrics()].
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     head(names(fn_schema()))
#'
#'     x <- fn_schema("adiv_boxplot")
#'     x$returns
#'     str(x$params[[10]])
#'
#'     fn_schema("rarefy", json = TRUE)
#'
fn_schema <- function (fn = NULL, json = FALSE) {

  validate_bool('json')

  if (is.null(fn)) {
    x <- lapply(setNames(nm = names(FN_SCHEMA)), schema_entry)
  } else {
    validate_schema_fn()
    x <- schema_entry(fn)
  }

  if (isTRUE(json)) return (schema_json(x))
  return (x)
}




#' Check arguments for an rbiom function against a dataset.
#'
#' Applies the checks the function itself would apply on entry, using its
#' [fn_schema()] entry, without running it. The dataset can be an rbiom
#' object or the summary returned by `fn_check_data()`, so that a program
#' that only keeps a summary can check requests before passing them on.
#'
#' Checks cover the type, choices, limits, and number of each value;
#' metadata field names and their classes; ranks; metric names (and
#' whether a phylogenetic metric has a tree); and the constraints in the
#' entry's `rules`. Checks that need computed data are left to the
#' function: names in `colors`, `shapes`, and `patterns` against the
#' groups, palette and color names, the parts of a `tracks` or `grid`
#' list, and taxa names.
#'
#' @param fn   The name of a function covered by [fn_schema()].
#'
#' @param args   A named list of the arguments to check, as they would be
#'        passed to `fn`. Arguments not given take their defaults.
#'        Arguments marked `r_only` in the schema can't be given.
#'
#' @param data   An [rbiom object][rbiom_objects], or the list returned by
#'        `fn_check_data()` (also after a round trip through JSON).
#'
#' @param biom   An [rbiom object][rbiom_objects], or any value accepted by
#'        [as_rbiom()].
#'
#' @return `fn_check()`: the arguments in `args`, converted as the
#'         function would convert them: abbreviations expanded, rank
#'         numbers replaced by rank names, layer letters by layer names,
#'         metric names by metric ids. A field given with a `"=="` or
#'         `"!="` prefix is returned without it, and added to `within` or
#'         `between`, which are then returned too. An invalid argument is
#'         an error, with the message the function would give.
#'
#'         `fn_check_data()`: a list with `fields` (one entry per
#'         metadata field, giving its `class`: `"numeric"` or
#'         `"categorical"`), `ranks` (the rank names), and `tree` (whether
#'         `biom` has a phylogenetic tree).
#'
#' @seealso [fn_schema()]
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     fn_check("adiv_boxplot", list(x = "body", adiv = "shan", layers = "xd"), hmp50)
#'
#'     data <- fn_check_data(hmp50)
#'     str(data)
#'
#'     fn_check("taxa_boxplot", list(rank = 2, taxa = 3), data)
#'
#'     try(fn_check("adiv_boxplot", list(ci = "zz"), data))
#'
fn_check <- function (fn, args = list(), data) {

  validate_schema_fn()

  if (!is.list(args) || is.data.frame(args))
    cli_abort("`args` must be a named list, not {.type {args}}.")

  if (length(args) > 0 && (is.null(names(args)) || !all(nzchar(names(args)))))
    cli_abort("Every element of `args` must be named.")

  specs <- FN_SCHEMA[[fn]]$params

  if (length(bad <- setdiff(names(args), names(specs))) > 0)
    cli_abort("{.fn {fn}} has no argument{?s} {.arg {bad}}.")

  r_only <- names(args)[vapply(specs[names(args)], function (s) isTRUE(s$r_only), logical(1))]
  if (length(r_only) > 0)
    cli_abort("{.arg {r_only}} can only be set by an R caller.")


  #________________________________________________________
  # The arguments as the function would see them: given ones, then
  # defaults (which may refer to other arguments, e.g. `stat.by = x`).
  #________________________________________________________
  fmls <- formals(get(fn, envir = asNamespace('rbiom'), mode = 'function'))
  env  <- rlang::new_environment(parent = rlang::ns_env('rbiom'))

  env$biom <- check_data_standin(data)

  for (i in setdiff(names(fmls), c('biom', '...'))) {
    if (hasName(args, i)) {
      env[[i]] <- args[[i]]
    } else if (isTRUE(specs[[i]]$r_only)) {
      next
    } else if (identical(fmls[[i]], quote(expr = ))) {
      cli_abort("{.fn {fn}} needs {.arg {i}}.")
    } else {
      do.call(delayedAssign, list(i, fmls[[i]], env, env))
    }
  }

  scratch <- schema_check(fn, env)


  #________________________________________________________
  # Return the given arguments, converted.
  #________________________________________________________
  out <- list()
  for (i in names(args)) {
    x <- strip_validated(get(i, envir = scratch, inherits = FALSE))
    if (identical(specs[[i]]$type, 'layers')) x <- unname(x)
    out[i] <- list(x)
  }

  for (i in c('within', 'between'))
    if (exists(i, envir = scratch, inherits = FALSE) && !eq(scratch[[i]], env[[i]]))
      out[i] <- list(scratch[[i]])

  return (out)
}


#' @rdname fn_check
#' @export
fn_check_data <- function (biom) {

  biom <- as_rbiom(biom)

  list(
    fields = lapply(biom$metadata, function (x) {
      list(class = if (is.numeric(x)) "numeric" else "categorical") }),
    ranks  = biom$ranks,
    tree   = !is.null(biom$tree) )
}




#________________________________________________________
# `fn` must name a covered function.
#________________________________________________________
validate_schema_fn <- function (env = parent.frame()) {

  fn <- get('fn', envir = env, inherits = FALSE)

  if (!is_scalar_character(fn) || is.na(fn))
    cli_abort("`fn` must be a function name, not {.type {fn}}.")

  if (!hasName(FN_SCHEMA, fn))
    cli_abort(c(
      'x' = "{.fn {fn}} is not described by {.fn fn_schema}.",
      'i' = "Covered functions: {.val {names(FN_SCHEMA)}}." ))

  return (invisible(NULL))
}


#________________________________________________________
# What the checks need from a dataset: `$metadata` (column names and
# classes), `$ranks`, and `$tree`. An rbiom object, or a stand-in built
# from fn_check_data().
#________________________________________________________
check_data_standin <- function (data) {

  if (inherits(data, 'rbiom')) return (data)

  if (!is.list(data) || !all(c('fields', 'ranks') %in% names(data)))
    cli_abort("`data` must be an rbiom object or a list from {.fn fn_check_data}.")

  fields <- as.list(data$fields)

  metadata <- lapply(fields, function (f) {
    if (identical(f[['class']], 'numeric')) numeric(0) else character(0) })

  list(
    metadata = structure(metadata, class = 'data.frame', row.names = integer(0)),
    ranks    = as.character(unlist(data$ranks)),
    tree     = if (isTRUE(as.logical(data$tree))) TRUE else NULL )
}


#________________________________________________________
# validate_bdiv() marks its result; callers don't need to see that.
#________________________________________________________
strip_validated <- function (x) {
  if (!is.null(attr(x, 'validated'))) attr(x, 'validated') <- NULL
  return (x)
}




#________________________________________________________
# The exported entry for one function.
#________________________________________________________
SCHEMA_FIELDS <- c(
  'type', 'required', 'r_only', 'deprecated', 'unused',
  'class', 'coerce', 'ops', 'extra', 'div', 'aes', 'by', 'forms', 'palettes', 'needs',
  'choices', 'choices_when', 'min', 'max', 'max_n', 'per_side', 'null',
  'default', 'default_expr', 'default_note', 'when', 'note', 'desc' )

schema_entry <- function (fn) {

  def   <- FN_SCHEMA[[fn]]
  fmls  <- formals(get(fn, envir = asNamespace('rbiom'), mode = 'function'))
  descs <- rd_arg_descs(fn)

  params <- lapply(names(def$params), function (nm) {

    spec <- def$params[[nm]]

    # Defaults, from formals(). (An argument without one holds the empty
    # symbol, which can't be assigned to a variable.)
    if (nm != '...') {
      if (identical(fmls[[nm]], quote(expr = ))) {
        spec$required <- TRUE
      } else if (is_literal(fml <- fmls[[nm]])) {
        spec['default'] <- list(eval(fml, baseenv()))
      } else {
        spec['default'] <- list(NULL)
        spec$default_expr <- paste(deparse(fml), collapse = " ")
      }
    }

    # Choices with display names.
    if (identical(spec$type, 'metric'))
      spec$choices <- metric_choices(spec$div)

    if (identical(spec$type, 'layers'))
      spec$choices <- unname(Map(
        f      = function (letter, name) list(letter = letter, name = name),
        letter = names(spec$choices),
        name   = unname(spec$choices) ))

    if (is.character(spec$choices))
      spec$choices <- unname(spec$choices)

    spec$desc <- descs[[nm]] %||% ""

    c(list(name = nm), spec[intersect(SCHEMA_FIELDS, names(spec))])
  })

  c(
    list(fn = fn, returns = def$returns),
    if (!is.null(def$rules)) list(rules = def$rules),
    list(params = params) )
}


#________________________________________________________
# A default that can be written as a value: a constant, or c() / unary
# minus of constants. Anything else (`x`, `n_cpus()`) is an expression.
#________________________________________________________
is_literal <- function (x) {
  if (is.null(x) || (is.atomic(x) && !is.symbol(x))) return (TRUE)
  if (!is.call(x)) return (FALSE)
  fun <- x[[1]]
  if (!is.symbol(fun) || !as.character(fun) %in% c('c', '-')) return (FALSE)
  if (as.character(fun) == '-' && length(x) != 2) return (FALSE)
  all(vapply(as.list(x)[-1], is_literal, logical(1)))
}


#________________________________________________________
# ecodive's metrics, as id / name / phylo / int_only.
#________________________________________________________
metric_choices <- function (div) {
  m <- ecodive::list_metrics(div)
  lapply(seq_len(nrow(m)), function (i) list(
    id       = m$id[[i]],
    name     = m$name[[i]],
    phylo    = isTRUE(m$phylo[[i]]),
    int_only = isTRUE(m$int_only[[i]]) ))
}


#________________________________________________________
# JSON with the list fields always written as arrays.
#________________________________________________________
SCHEMA_ARRAYS <- c(
  'choices', 'forms', 'ops', 'extra', 'palettes', 'has',
  'is', 'is_not', 'has_any', 'one_of', 'any_of', 'requires', 'forbids' )

schema_json <- function (x) {

  arrays <- function (x) {
    if (!is.list(x)) return (x)
    for (i in seq_along(x)) {
      xi <- x[[i]]
      if (is.null(xi)) next
      if (!is.null(names(x)) && names(x)[[i]] %in% SCHEMA_ARRAYS && is.atomic(xi)) {
        x[[i]] <- I(xi)
      } else if (is.list(xi)) {
        x[[i]] <- arrays(xi)
      }
    }
    x
  }

  jsonlite::toJSON(arrays(x), auto_unbox = TRUE, null = "null", digits = NA)
}


#________________________________________________________
# One-line descriptions of a function's arguments: the first sentence of
# each argument's text in the help page.
#________________________________________________________
SCHEMA_CACHE <- new.env(parent = emptyenv())

rd_arg_descs <- function (fn) {

  SCHEMA_CACHE$descs <- SCHEMA_CACHE$descs %||% list()
  if (!is.null(SCHEMA_CACHE$descs[[fn]])) return (SCHEMA_CACHE$descs[[fn]])

  rd <- rd_for(fn)
  if (is.null(rd)) return (list())

  descs <- list()
  for (sec in rd) {
    if (!identical(attr(sec, 'Rd_tag'), '\\arguments')) next
    for (item in sec) {
      if (!identical(attr(item, 'Rd_tag'), '\\item') || length(item) < 2) next
      nms  <- trimws(strsplit(rd_text(item[[1]]), ",")[[1]])
      text <- first_sentence(rd_text(item[[2]]))
      for (nm in nms) descs[[nm]] <- text
    }
  }

  SCHEMA_CACHE$descs[[fn]] <- descs
  return (descs)
}


#________________________________________________________
# The parsed help page that documents `fn`.
#________________________________________________________
rd_for <- function (fn) {

  if (is.null(SCHEMA_CACHE$rd_db)) {

    # Under pkgload::load_all(), system.file() points into the source
    # tree's `inst` directory; an installed package has no `man`.
    pkg <- system.file(package = 'rbiom')
    src <- Filter(
      f = function (d) all(file.exists(file.path(d, c('DESCRIPTION', 'man')))),
      x = c(pkg, dirname(pkg)) )

    SCHEMA_CACHE$rd_db <- tryCatch(
      expr  = if (length(src) > 0) tools::Rd_db(dir = src[[1]]) else tools::Rd_db('rbiom'),
      error = function (e) list() )
  }

  for (rd in SCHEMA_CACHE$rd_db) {
    aliases <- unlist(lapply(rd, function (x) {
      if (identical(attr(x, 'Rd_tag'), '\\alias')) rd_text(x) }))
    if (fn %in% trimws(aliases)) return (rd)
  }

  return (NULL)
}


#________________________________________________________
# Plain text of an Rd fragment, with code in backticks.
#________________________________________________________
rd_text <- function (x) {

  tag <- attr(x, 'Rd_tag') %||% ""

  if (tag %in% c('COMMENT', '\\dontrun', '\\dontshow', '\\testonly')) return ("")
  if (tag == '\\dots' || tag == '\\ldots') return ("...")
  if (tag == '\\R') return ("R")
  if (tag == '\\cr') return (" ")
  if (tag == '\\href' && length(x) >= 2) return (rd_text(x[[2]]))

  txt <- if (is.list(x)) paste(vapply(x, rd_text, character(1)), collapse = "") else paste(x, collapse = "")

  if (tag %in% c('\\code', '\\verb', '\\env', '\\file', '\\samp', '\\option'))
    txt <- paste0("`", txt, "`")

  if (tag == '\\sQuote') txt <- paste0("'", txt, "'")
  if (tag == '\\dQuote') txt <- paste0('"', txt, '"')

  return (txt)
}


#________________________________________________________
# The first sentence of a paragraph, whitespace squished. A very short
# first sentence (e.g. "Logical.") takes the next one with it.
#________________________________________________________
first_sentence <- function (txt) {
  txt  <- trimws(gsub("\\s+", " ", txt))
  ends <- gregexpr("(?<!\\be\\.g|\\bi\\.e|\\bvs|\\betc)\\.(?=\\s|$)", txt, perl = TRUE)[[1]]
  ends <- ends[ends > 0]
  end  <- ends[ends >= 20][1]
  if (is.na(end)) end <- nchar(txt)
  return (substr(txt, 1, end))
}
