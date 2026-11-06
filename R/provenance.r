#________________________________________________________
# Provenance tracking.
#
# Two mechanisms work together so that every rbiom result can be regenerated
# from code:
#
# 1. rbiom objects carry a *history*: the list of calls (constructor, then
#    one entry per transformation) that produced the object's current state.
#    `biom$code` prints it. Transformation functions call `prov_quiet()`
#    around their internal edits and then `prov_log()` once with the
#    user-level call. Direct assignments (`biom$metadata <- x`) are logged
#    by the active-binding setters as entries that are flagged as not
#    reproducible from code.
#
# 2. Outputs (tables, stats, plots) carry `cmd` (the rbiom call that made
#    them), `history` (a snapshot of the input object's code), and for stats
#    and plots `code` (the statistics/ggplot2 code that rebuilds them from
#    `data`). `assemble_code()` joins these into the `$code` shown to users.
#    Because provenance depends on the object's history and not only on its
#    contents, these attributes are attached *after* any cache lookup - see
#    `prov_attach()`.
#
# 3. `P('pkg::fn')` wraps functions so that their calls are recorded as
#    display strings for the generated ggplot2 code.
#________________________________________________________

PROV <- new.env()
PROV$depth <- 0L



#________________________________________________________
# Run `expr` without the setters logging history entries.
#________________________________________________________
prov_quiet <- function (expr) {
  PROV$depth <- PROV$depth + 1L
  on.exit(PROV$depth <- PROV$depth - 1L, add = TRUE)
  expr
}


#________________________________________________________
# Read / replace an rbiom object's history entries.
#________________________________________________________
prov_get <- function (biom) {
  biom$.__enclos_env__$private$.history
}

prov_set <- function (biom, entries) {
  biom$.__enclos_env__$private$.history <- entries
  invisible(biom)
}


#________________________________________________________
# Append a history entry to `biom`.
#
# call         - A call in which the symbol `biom` denotes the object's
#                previous state, e.g. quote(subset(biom, Age > 30)). May
#                also be a string that is shown verbatim.
# reproducible - FALSE when the call can't be re-run as written (because
#                it references in-memory data). Such lines are shown as
#                comments in `$code`.
# note         - Optional trailing comment.
# prefix       - Whether to print the entry as `biom <- <call>` (TRUE) or
#                as-is (FALSE, for assignments like `biom$id <- "x"`).
# inputs       - For biom_merge(): a list of the input objects' histories.
#________________________________________________________
prov_log <- function (biom, call, reproducible = TRUE, note = NULL, prefix = TRUE, inputs = NULL) {

  if (is.character(call)) call <- paste(call, collapse = "\n")

  entry <- list(
    call         = call,
    reproducible = isTRUE(reproducible),
    note         = note,
    prefix       = isTRUE(prefix),
    inputs       = inputs,
    hash         = biom$hash )

  prov_set(biom, c(prov_get(biom), list(entry)))
}


#________________________________________________________
# Replace the history with a single constructor entry.
# `pc` is the result of `prov_call()`.
#________________________________________________________
prov_init <- function (biom, pc) {
  prov_set(biom, list())
  prov_log(biom, call = pc$call, reproducible = pc$reproducible)
}


#________________________________________________________
# Called by the active-binding setters. Logs a (usually non-reproducible)
# assignment unless we're inside an rbiom transformation.
#________________________________________________________
prov_setter <- function (biom, field, value) {

  if (PROV$depth > 0L) return (invisible(NULL))

  pv   <- prov_value(value)
  lhs  <- call('$', as.name('biom'), as.name(field))
  call <- as.call(list(as.name('<-'), lhs, pv$value))

  prov_log(biom, call = call, reproducible = pv$reproducible, prefix = FALSE)
}


#________________________________________________________
# A short description of an object, for history entries that can't
# embed the object itself.
#________________________________________________________
obj_repr <- function (x) {

  if (is.null(x))                 return ("NULL")
  if (inherits(x, 'rbiom'))       return ("<rbiom object>")
  if (inherits(x, 'phylo'))       return (sprintf("<phylo; %d tips>", length(x$tip.label)))
  if (is.data.frame(x))           return (sprintf("<%s; %d x %d>", class(x)[[1]], nrow(x), ncol(x)))
  if (is.matrix(x))               return (sprintf("<matrix; %d x %d>", nrow(x), ncol(x)))
  if (inherits(x, 'Matrix'))      return (sprintf("<%s; %d x %d>", class(x)[[1]], nrow(x), ncol(x)))
  if (is.function(x))             return ("<function>")
  if (is.atomic(x) && length(x) == 1 && !is.character(x)) return (deparse(x))
  if (is.character(x) && length(x) == 1) return (sprintf("<string; %d chars>", nchar(x)))
  if (is.atomic(x))               return (sprintf("<%s [%d]>", class(x)[[1]], length(x)))
  if (is.list(x))                 return (sprintf("<list [%d]>", length(x)))

  return (sprintf("<%s>", class(x)[[1]]))
}


#________________________________________________________
# Represent a value for embedding in a logged call. Small atomic values
# are embedded as-is; anything else becomes a `<description>` placeholder
# and marks the entry as not reproducible.
#________________________________________________________
prov_value <- function (x, max_len = 50L, max_chars = 500L) {

  if (is.null(x)) return (list(value = NULL, reproducible = TRUE))

  small <- is.atomic(x) &&
    length(x) <= max_len &&
    is.null(dim(x)) &&
    (!is.character(x) || all(nchar(x) <= max_chars, na.rm = TRUE)) &&
    (is.null(attributes(x)) || identical(names(attributes(x)), 'names'))

  if (small) return (list(value = x, reproducible = TRUE))

  list(value = as.name(obj_repr(x)), reproducible = FALSE)
}


#________________________________________________________
# Build a call `fn(args...)` with each argument passed through
# `prov_value()`. Returns list(call, reproducible).
#________________________________________________________
prov_call <- function (fn, args = list()) {

  ok   <- TRUE
  vals <- lapply(args, function (a) {
    if (is.language(a)) return (a)
    pv <- prov_value(a)
    if (!pv$reproducible) ok <<- FALSE
    pv$value
  })

  list(
    call         = as.call(c(as.name(fn), vals)),
    reproducible = ok )
}


#________________________________________________________
# Inline the values of free variables in a non-standard-evaluation
# expression (e.g. the `subset` argument), so that `Age > cutoff` is
# recorded as `Age > 30`. Column names (`exclude`) and functions are left
# alone. Variables whose values can't be embedded make the result
# non-reproducible.
#________________________________________________________
prov_inline <- function (expr, env, exclude = character()) {

  ok      <- TRUE
  exclude <- c(exclude, '.data', '.env', 'T', 'F')

  walk <- function (e) {

    if (is.symbol(e)) {
      nm <- as.character(e)
      if (!nzchar(nm) || nm %in% exclude)       return (e)
      if (!exists(nm, envir = env))             return (e)
      val <- get(nm, envir = env)
      if (is.function(val))                     return (e)
      pv <- prov_value(val)
      if (!pv$reproducible) { ok <<- FALSE;     return (e) }
      return (pv$value)
    }

    if (is.call(e)) {

      head <- e[[1]]
      if (is.symbol(head) && as.character(head) %in% c('function', '::', ':::', '~')) return (e)

      # For `x$name` and `x@slot`, only the object can be a variable.
      last <- if (is.symbol(head) && as.character(head) %in% c('$', '@')) 2L else length(e)

      for (i in seq_len(last)[-1])
        if (!identical(e[[i]], quote(expr = )))
          e[[i]] <- walk(e[[i]])

      return (e)
    }

    return (e)
  }

  list(expr = walk(expr), reproducible = ok)
}


#________________________________________________________
# R literals for a double vector that parse back to exactly the same values
# on every platform.
#
# deparse()'s "digits17" option is not enough: R's parser accumulates the
# digits of a decimal literal in a long double, which on some platforms
# (macOS on Apple silicon) is no wider than a double, so a 17-digit literal
# such as 15.299999999999999 can come back one unit in the last place off.
# Up to 15 significant digits with a power-of-ten scale within +/-22 always
# parse exactly (the digits fit in 53 bits and the power of ten is exact),
# so that form is used whenever it round-trips. Anything else - typically
# floating-point noise such as 9 * 1.7 - is written in C99 hexadecimal
# notation (0x1.e99999999999ap+3), which R parses exactly everywhere.
#________________________________________________________
num_literal <- function (x) {

  x   <- as.double(x)
  out <- sprintf("%.15g", x) # also gives "NA", "NaN", "Inf", "-Inf"
  fin <- is.finite(x)

  if (any(fin)) {
    s     <- sub("^-", "", out[fin])
    mant  <- sub("e.*$", "", s)
    expo  <- ifelse(grepl("e", s, fixed = TRUE), sub("^.*e", "", s), "0")
    scale <- as.integer(expo) - nchar(sub("^[^.]*\\.?", "", mant))
    exact <- abs(scale) <= 22L & as.numeric(out[fin]) == x[fin]
    out[fin][!exact] <- sprintf("%a", x[fin][!exact])

    # R doesn't read subnormal hexadecimal literals back exactly.
    bad <- as.numeric(out[fin]) != x[fin]
    out[fin][bad] <- sprintf("%.17g", x[fin][bad])
  }

  if (length(x) > 0 && all(is.na(x) & !is.nan(x))) out[] <- "NA_real_"

  return (out)
}


# TRUE when deparse()'s default (15 significant digits) can't represent `x`
# exactly on every platform.
needs_num_literal <- function (x) {
  is.double(x) && length(x) > 0 && any(grepl("0x", num_literal(x), fixed = TRUE))
}


#________________________________________________________
# deparse() for generated code, with every double constant written by
# num_literal(). Double vectors that need hexadecimal digits are swapped for
# placeholder symbols, the expression is deparsed as usual (so names and
# strings are still quoted and escaped by deparse()), and the placeholders
# are then replaced with the literals. The replacement text is built from
# numbers only.
#________________________________________________________
deparse_exact <- function (expr, width.cutoff = 60L, control = c("keepNA", "keepInteger", "niceNames", "showAttributes")) {

  lits <- character()
  pfx  <- "RBIOMNUM"
  while (any(grepl(pfx, deparse(expr, control = NULL), fixed = TRUE)))
    pfx <- paste0(pfx, "X")

  literal <- function (x) {
    lit <- num_literal(x)
    if (length(x) == 1 && is.null(names(x))) return (lit)
    call <- as.call(c(as.name("c"), lapply(seq_along(x), function (i) as.name(sprintf("%s%dV_", pfx, i)))))
    if (!is.null(names(x))) names(call) <- c("", names(x))
    txt <- paste(deparse(call, width.cutoff = width.cutoff, control = control), collapse = "\n")
    for (i in rev(seq_along(x)))
      txt <- gsub(sprintf("%s%dV_", pfx, i), lit[[i]], txt, fixed = TRUE)
    txt
  }

  walk <- function (e) {

    if (needs_num_literal(e) && is.null(dim(e)) &&
        (is.null(attributes(e)) || identical(names(attributes(e)), 'names'))) {
      lits[[length(lits) + 1L]] <<- literal(e)
      return (as.name(sprintf("%s%dE_", pfx, length(lits))))
    }

    if (is.call(e)) {
      for (i in seq_along(e))
        if (!identical(e[[i]], quote(expr = )) && !is.null(e[[i]]))
          e[[i]] <- walk(e[[i]])
    }

    return (e)
  }

  if (!is.call(expr) && !needs_num_literal(expr))
    return (deparse(expr, width.cutoff = width.cutoff, control = control))

  txt <- deparse(walk(expr), width.cutoff = width.cutoff, control = control)
  for (k in rev(seq_along(lits)))
    txt <- gsub(sprintf("%s%dE_", pfx, k), lits[[k]], txt, fixed = TRUE)

  return (txt)
}


#________________________________________________________
# Format history entries as lines of R code.
#________________________________________________________
prov_entry_lines <- function (entries, var = "biom") {

  lines <- character()

  for (e in entries) {

    for (k in seq_along(e$inputs))
      lines <- c(lines, prov_entry_lines(e$inputs[[k]], paste0(var, k)), "")

    if (is.character(e$call)) {
      txt <- e$call

    } else {
      # deparse() quotes and escapes embedded values, so data-originating
      # strings (sample names, metadata levels, paths) can't inject code.
      # Placeholders such as `<tbl_df; 50 x 5>` keep their backticks; such
      # entries are always emitted as comments.
      call <- do.call(substitute, list(e$call, list(biom = as.name(var))))
      txt  <- paste(deparse_exact(call, width.cutoff = 120L), collapse = "\n  ")
    }

    line <- if (isTRUE(e$prefix)) sprintf("%s <- %s", var, txt) else txt

    if (!isTRUE(e$reproducible)) {
      line <- paste0("# ", gsub("\n", "\n# ", line), "   # not reproducible from code")
    } else if (!is.null(e$note)) {
      line <- paste0(line, "   # ", e$note)
    }

    lines <- c(lines, line)
  }

  return (lines)
}


#________________________________________________________
# The R code that reproduces `biom`, as a character vector of lines.
#________________________________________________________
prov_lines <- function (biom, var = "biom") {

  entries <- prov_get(biom)
  lines   <- prov_entry_lines(entries, var)

  if (length(entries) > 0 && !identical(entries[[length(entries)]]$hash, biom$hash))
    lines <- c(lines, sprintf("# %s was subsequently modified in a way that was not recorded.", var))

  return (lines)
}


#________________________________________________________
# Attach provenance to a result. Call this after any cache lookup.
#________________________________________________________
prov_attach <- function (res, cmd = NULL, biom = NULL, data = NULL) {

  if (is.null(res)) return (NULL)

  if (!is.null(cmd)) attr(res, 'cmd') <- as.character(cmd)

  if (inherits(biom, 'rbiom')) {
    attr(res, 'history') <- prov_lines(biom)

  } else if (!is.null(data)) {
    attr(res, 'history')  <- attr(data, 'history', exact = TRUE)
    attr(res, 'data_cmd') <- attr(data, 'cmd',     exact = TRUE)
  }

  return (res)
}


#________________________________________________________
# Compose the full `$code` for a table, stats, or plot object:
#
#   library(rbiom) [+ library(ggplot2)]
#   <history of the rbiom object>
#   data <- <cmd that produced the data>
#   <stats / ggplot2 code>
#________________________________________________________
assemble_code <- function (obj) {

  cmd      <- attr(obj, 'cmd',      exact = TRUE)
  code     <- attr(obj, 'code',     exact = TRUE)
  history  <- attr(obj, 'history',  exact = TRUE)
  data_cmd <- attr(obj, 'data_cmd', exact = TRUE)

  data <- tryCatch(obj[['data']], error = function (e) NULL)
  if (!is.data.frame(data)) data <- NULL

  if (is.null(history)  && !is.null(data)) history  <- attr(data, 'history', exact = TRUE)
  if (is.null(data_cmd) && !is.null(data)) data_cmd <- attr(data, 'cmd',     exact = TRUE)

  if (is.null(data_cmd) && is.null(code)) data_cmd <- cmd

  is_plot <- inherits(obj, c('ggplot', 'patchwork', 'gg'))

  lines <- c(
    sprintf("# Generated by rbiom %s", packageVersion('rbiom')),
    "library(rbiom)",
    if (is_plot) "library(ggplot2)",
    "" )

  if (length(history) > 0)  lines <- c(lines, history, "")
  if (!is.null(data_cmd))   lines <- c(lines, paste("data <-", data_cmd), "")
  if (!is.null(code))       lines <- c(lines, code)

  lines <- sub("\\s+$", "", lines)
  while (length(lines) > 0 && !nzchar(lines[[length(lines)]])) lines <- lines[-length(lines)]

  add_class(paste(lines, collapse = "\n"), 'rbiom_code')
}




#________________________________________________________
# Provenance-tracked wrappers for functions used to build ggplot2 code.
#________________________________________________________

FUNCS <- new.env()

P <- function (pkgfn) {

  if (env_has(FUNCS, pkgfn))
    return (FUNCS[[pkgfn]])


  parts <- strsplit(pkgfn, '::', fixed = TRUE)[[1]]
  pkg   <- parts[[1]]
  fn    <- parts[[2]]

  require_package(pkg, reason = paste0('to use `', fn, '()`'))

  fun <- getFromNamespace(x = fn, ns = pkg)
  fn  <- ifelse(pkg %in% c('base', 'ggplot2', 'grid', 'rbiom'), fn, pkgfn)


  newfun <- function (..., .indent = 0) {

    dots <- list(...)
    res  <- fun(...)

    argstr  <- as.args(dots, fun = fun, indent = .indent)
    display <- sprintf(fmt = "%s(%s)", fn, argstr)

    attr(res, 'display') <- display

    return (res)
  }

  attr(newfun, 'display')    <- fn
  attr(newfun, 'formalArgs') <- formalArgs(fun)
  attr(newfun, 'formals')    <- formals(fun)

  FUNCS[[pkgfn]] <- newfun

  return (newfun)
}
