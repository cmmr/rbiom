


#____________________________________________________________________
# Ensure something's a factor with minimal levels.
#____________________________________________________________________
refactor <- function (values) {
  values <- as.factor(values)
  factor(values, intersect(levels(values), unique(values)))
}


#____________________________________________________________________
# Assign var = value only when var doesn't already exist.
#____________________________________________________________________
default <- function (var, value, env = parent.frame()) {
  if (!hasName(env, var)) assign(var, value, env)
  return (invisible(NULL))
}


#____________________________________________________________________
# A string representation of the function call.
#____________________________________________________________________
current_cmd <- function (fn, indent = 0) {
  
  fun  <- get(fn)
  env  <- rlang::caller_env(2)
  args <- rlang::caller_call() %>% 
    rlang::call_match(fun) %>% 
    rlang::call_args() %>% 
    lapply(eval, envir = env)
  
  glue("{fn}({as.args(args, indent, fun)})")
}


#____________________________________________________________________
# Collect and evaluate all the arguments.
#____________________________________________________________________
slurp_env <- function (..., .dots = FALSE) {
  
  env <- rlang::caller_env()
  
  if (isTRUE(.dots)) { env$.dots <- lapply(list(...), eval, envir = env)
  } else             { rlang::env_bind(env, ...) }
  
  rlang::env_unbind(env, '...')
  params <- lapply(as.list(env, all.names = TRUE), eval, envir = env)
  rlang::env_unbind(env, names(env))
  
  return (params)
}



#____________________________________________________________________
# Create an environment of evaluated variables from caller's env.
#____________________________________________________________________
eval_envir <- function (env, ...) {
  
  dots <- list(...)
  
  
  #____________________________________________________________________
  # Evaluate all arguments into a clean env.
  #____________________________________________________________________
  params <- lapply(as.list(env), eval) %>%
    rlang::new_environment(parent = rlang::ns_env('rbiom'))
  
  
  #____________________________________________________________________
  # Check if obviously wrong parameters are going into dots.
  #____________________________________________________________________
  if (length(dots) > 0) {
    
    invalid <- paste(collapse = ", ", c(
      names(dots)[startsWith(names(dots), '.')],
      intersect(
        x = names(dots),
        y = c(
          'x', 'y', 'color.by', 'shape.by', 'facet.by', 
          'pattern.by', 'label.by', 'order.by', 'stat.by', 'limit.by',
          'adiv', 'bdiv', 'taxa', 'rank', 'weighted', 
          'within', 'between', 'tree', 'params',
          'regr', 'resp' ))))
    
    if (nzchar(invalid)) {
      fn <- attr(rlang::caller_fn(), 'fn') %||% "function"
      cli_abort(c(x = "{fn} does not accept parameter(s): {invalid}"))
    }
    
  }
  
  
  params$.dots <- lapply(dots, eval)
  
  
  return (params)
}


#____________________________________________________________________
# Check if two objects are identical, ignoring attr(,'hash').
#____________________________________________________________________
eq <- function (x, y) {
  if (hasName(attributes(x), 'hash')) attr(x, 'hash') <- NULL
  if (hasName(attributes(y), 'hash')) attr(y, 'hash') <- NULL
  identical(x, y)
}




#____________________________________________________________________
# Convert to tibble and add 'rbiom_tbl' class.
#____________________________________________________________________
as_rbiom_tbl <- function (df) {
  add_class(as_tibble(df), 'rbiom_tbl')
}



#____________________________________________________________________
# Explicitly define the code to be displayed in cmd
#____________________________________________________________________
as.cmd <- function (expr, env=NULL) {
  if (is_null(env)) {
    expr <- substitute(expr)
    cmd  <- capture.output(expr)
  } else {
    expr <- do.call(substitute, list(expr=substitute(expr), env=env))
    cmd  <- capture.output(expr)
  }
  
  if (length(cmd) > 1)
    cmd <- paste(trimws(cmd), collapse = " ")
  
  structure(eval(expr), 'display' = cmd)
}



#____________________________________________________________________
# Convert a list of arguments to character strings.
# When indent > 0, produces a multi-line string.
# When fun is a function, puts the arguments in the expected order.
#____________________________________________________________________
as.args <- function (args = parent.frame(), indent = 0, fun = NULL) {
  
  if (is.environment(args)) args <- as.list(args)
  stopifnot(is_list(args))
  stopifnot(is_scalar_integerish(indent) && !is_na(indent))
  
  
  # Discard arguments with `display = FALSE` attribute.
  for (i in rev(seq_along(args)))
    if (isFALSE(attr(args[[i]], 'display', exact = TRUE)))
      args[[i]] <- NULL
  
  
  # Re-arrange parameter order; omit names where implicitly known; ignore if default.
  if (is.function(fun) && !is_null(names(args))) {
    f_args <- formals(fun)
    
    for (i in names(args))
      if (hasName(f_args, i) && eq(args[[i]], f_args[[i]]))
        args[[i]] <- NULL
    
    f_args <- names(f_args)
    args   <- args[c(intersect(f_args, names(args)), sort(setdiff(names(args), f_args)))]
    if (!'...' %in% f_args) args <- args[intersect(f_args, names(args))]
    if (isTRUE(indent == 0))
      for (i in seq_along(args))
        if (eq(names(args)[[i]], f_args[[i]])) names(args)[i] <- "" else break
  }
  
  # Right-pad parameter names.
  fmt <- "%s = %s"
  if (isTRUE(indent > 0 && length(args) > 1 && !is_null(names(args))))
    fmt <- paste(collapse="", rep(" ", indent)) %>%
      paste0("\n", ., "%-", max(nchar(coan(names(args)))), "s = %s")
  
  # Convert arguments to `eval`-able string representations
  strs <- c()
  for (i in seq_along(args)) {
    
    key <- if (is_null(names(args))) '' else names(args)[[i]]
    val <- args[[i]]
    
    display <- attr(val, 'display', exact = TRUE)
    
    val <- if (is_null(val))          { "NULL"
    } else if (!is_null(display))     { display
    } else if (is_character(val))     { double_quote(val) 
    } else if (is_logical(val))       { as.character(val) %>% setNames(names(val))
    } else if (is.numeric(val))       { as.character(val) %>% setNames(names(val))
    } else if (inherits(val, 'quosures'))   { as.character(val)
    } else if (inherits(val, 'rbiom'))      { "biom"
    } else if (is.data.frame(val))    { "data"
    } else if (inherits(val, 'formula'))    { format(val)
    } else if (is.function(val))      { fun_toString(val)
    } else if (inherits(val, 'uneval'))     { aes_toString(val)
    } else if (is.factor(val))        { as.character(val) %>% setNames(names(val))
    } else if (is_list(val))          { paste0("list(", as.args(val), ")")
    } else                            { capture.output(val) }
    
    
    if (isTRUE(is_character(val) && !is_null(names(val))))
      val <- paste(single_quote(names(val)), "=", unname(val))
    
    if (length(val) > 1)
      val <- paste0("c(", paste(collapse=", ", val), ")")
    
    if (nzchar(key)) {
      val <- sprintf(fmt = fmt, coan(key), val)
    }
    
    strs <- c(strs, val)
  }
  
  strs <- paste(strs, collapse = ", ")
  
  if (isTRUE(indent > 0 && length(args) > 1))
    strs <- paste0(strs, " ")
  
  return (strs)
}


#____________________________________________________________________
# Find a function's name
#____________________________________________________________________
fun_toString <- function (x) {
  
  pkg <- environment(x)[['.packageName']]
  if (!is_null(pkg))
    for (fn in getNamespaceExports(pkg))
      if (eq(x, do.call(`::`, list(pkg, fn)))) {
        if (!pkg %in% getOption("defaultPackages"))
          fn <- sprintf("%s::%s", pkg, fn)
        return (fn)
      }
  
  if (!is.null(srcref <- attr(x, 'srcref',  exact = TRUE)))
    return (capture.output(srcref))
  
  chr <- paste0(capture.output(x), collapse = "")
  if (nchar(chr) < 50) return (chr)
  
  return ("_Custom_Function")
}


#____________________________________________________________________
# Convert an aes object to a string
#____________________________________________________________________
aes_toString <- function (x) {
  
  # Consistently order the aes parameters
  keys <- names(x)
  sort <- c(
    "x", "xend", "xmin", "xmax", 
    "y", "yend", "ymin", "ymax", 
    "colour", "fill", "shape", "group", "label" )
  keys <- c(intersect(sort, keys), setdiff(keys, sort))
  
  results <- c()
  for (key in keys) {
    
    val <- x[[key]]
    val <- if (inherits(val, 'quosure')) { capture.output(rlang::quo_get_expr(val))
    } else if (inherits(val, 'formula')) { capture.output(as.name(all.vars(val)))
    } else if (is.logical(val))    { as.character(val)
    } else if (is.numeric(val))    { as.character(val)
    } else                         { double_quote(val) }
    
    key %<>% sub(pattern = "colour", replacement = "color")
    key <- capture.output(as.name(key))
    
    results %<>% c(sprintf("%s = %s", key, val))
  }
  
  return (sprintf("aes(%s)", paste(collapse = ", ", results)))
}


#____________________________________________________________________
# A command as a string. Assumes dist=dist, ord=ord, etc.
# Example: fmt_cmd("dm <- %s", bdiv_distmat, biom, weighted, tree)
#____________________________________________________________________
fmt_cmd <- function (.fmt, .fun, ...) {
  x    <- all.vars(match.call())[-c(1,2)]
  fn   <- deparse(substitute(.fun))
  env  <- parent.frame()
  args <- sapply(x, simplify = FALSE, function (i) {
    if (endsWith(i, "_"))
      return (aa(NA, display = substr(i, 1, nchar(i) - 1)))
    get(i, envir = env, inherits = TRUE)
  })
  names(args) %<>% sub("_$", "", .)
  sprintf(.fmt, sprintf("%s(%s)", fn, as.args(args, fun = .fun)))
}


#____________________________________________________________________
# Add attributes to an object.
#____________________________________________________________________
aa <- function (obj, ...) {
  dots <- list(...)
  if (!is.null(obj))
    for (k in names(dots))
      attr(obj, k) <- dots[[k]]
  return (obj)
}


#____________________________________________________________________
# Drop/Keep particular columns from a data.frame
#____________________________________________________________________
drop_cols <- function (df, ...) {
  cols <- unlist(list(...))
  if (is_null(cols) || is_null(df)) return (df)
  df[,setdiff(colnames(df), cols),drop=FALSE]
}
keep_cols <- function (df, ...) {
  cols <- unique(unlist(list(...)))
  if (is_null(cols) || is_null(df)) return (NULL)
  df[,intersect(cols, colnames(df)),drop=FALSE]
}
rename_cols <- function(df, ...) {
  
  vals <- list(...)
  if (length(vals) == 1 && is.list(vals[[1]]))
    vals <- as.list(vals[[1]])
  
  x <- attr(df, 'names', exact = TRUE)
  for (i in names(vals))
    if (i %in% x)
      x[which(x == i)] <- vals[[i]]
  
  attr(df, 'names') <- x
  return (df)
}


#____________________________________________________________________
# Decorate column names for plyr
#____________________________________________________________________
ply_cols <- function (cols) {
  
  if (inherits(cols, 'quoted')) return (cols)
  
  structure(lapply(cols, as.name), class = "quoted")
  
  # vars <- NULL
  # for (col in cols)
  #   vars <- c(plyr::as.quoted(as.name(col)), vars)
  # 
  # return (vars)
}





#____________________________________________________________________
# Log scale labels
#____________________________________________________________________
loglabels <- function (values) {
  
  hi <- ceiling(log10(max(force(values))))
  list(
    'breaks'       = as.cmd(10 ** (0:hi),                    env = list(hi = hi)),
    'minor_breaks' = as.cmd(as.vector(2:9 %o% 10 ** (0:hi)), env = list(hi = hi - 1)),
    'labels'       = label_number(scale_cut = cut_si("")) )
}

# siunit <- function (x) {
#   sapply(as.numeric(x), function (n) {
#     if (is.na(n))       return ("")
#     if (abs(n) < 10^3 ) return (as.character(n))
#     if (abs(n) < 10^6 ) return (sprintf("%s k", round(n / 10^3,  1)))
#     if (abs(n) < 10^9 ) return (sprintf("%s M", round(n / 10^6,  1)))
#     if (abs(n) < 10^12) return (sprintf("%s G", round(n / 10^9,  1)))
#     if (abs(n) < 10^15) return (sprintf("%s T", round(n / 10^12, 1)))
#     if (abs(n) < 10^18) return (sprintf("%s P", round(n / 10^15, 1)))
#     if (abs(n) < 10^21) return (sprintf("%s E", round(n / 10^18, 1)))
#     return (as.character(scientific(n)))
#   })
# }

#____________________________________________________________________
# Easily create list(x = ".x", label = ".label")
#____________________________________________________________________
.qw <- function (...) {
  x <- all.vars(match.call())
  as.list(setNames(paste0(".", x), x))
}


#____________________________________________________________________
# Capture Output As Name
#____________________________________________________________________
coan <- function (nm) {
  if (is.null(nm)) return ('')
  sapply(nm, USE.NAMES = FALSE, function (nm) capture.output(as.name(nm)))
}


#____________________________________________________________________
# Convert vector to a string, limited by width.
#____________________________________________________________________
vw <- function (x, width = min(80, floor(getOption("width") * 0.75))) {
  
  if (length(x) == 0) return ("<none>")
  if (length(x) == 1) return (x)
  if (length(x) == 2) return (paste(x[[1]], "and", x[[2]]))
  
  x[length(x)] <- paste0("and ", x[length(x)])
  x[-1]        <- paste0(", ",   x[-1])
  
  if (sum(nchar(x)) <= width)
    return (paste(x, collapse = ''))
  
  n <- sum(cumsum(nchar(x)) < width - nchar(", ..."))
  
  if (n <= 3)
    return (paste(c(head(x, max(1, n)), ", ..."), collapse = ''))
  
  n <- sum(cumsum(nchar(x)) < width - nchar(", ...") - nchar(tail(x, 1)))
  if (n <= 6)
    return (paste(c(head(x, n), ", ...", tail(x, 1)), collapse = ''))
  
  n <- sum(cumsum(nchar(x)) < width - nchar(", ...") - sum(nchar(tail(x, 2))))
  return (paste(c(head(x, n), ", ...", tail(x, 2)), collapse = ''))
  
}


#____________________________________________________________________
# Prepend a new class to an object's class list.
#____________________________________________________________________
add_class <- function (obj, cls) {
  if (!is.null(obj))
    class(obj) <- c(cls, setdiff(class(obj), cls))
  return (obj)
}





#____________________________________________________________________
# Shorthand for replacing NULLs with something else.
#____________________________________________________________________
if.null <- function (x, replacement) {
  if (is_null(x)) replacement else x
}


#____________________________________________________________________
# Like ifelse() but test is length 1 vector
#____________________________________________________________________
bool_switch <- function (test, yes, no) {
  res <- if (isTRUE(test)) yes else no
  return (res)
}


#____________________________________________________________________
# Find which parameters can be passed to a function.
#____________________________________________________________________
fun_params <- function (fun, params) {
  
  if (is.environment(params))
    params <- as.list(params, all.names = TRUE)
  
  stopifnot(is.function(fun))
  stopifnot(is.list(params))
  
  if (is.list(params$.dots))
    params <- c(params[names(params) != ".dots"], params$.dots)
  
  if (!"..." %in% formalArgs(fun))
    params <- params[intersect(formalArgs(fun), names(params))]
  
  return (params)
}




#____________________________________________________________________
# Remove unused factor levels and whitespace from a tibble.
#____________________________________________________________________
relevel <- function (tbl) {
  
  for (i in seq_len(ncol(tbl))) {
    
    # Enforce unnamed vectors.
    if (!is.null(names(tbl[[i]]))) tbl[[i]] %<>% unname()
    
    if (is.factor(tbl[[i]])) {
      
      lvls <- levels(tbl[[i]])
      lvls[which(!(lvls %in% as.character(tbl[[i]])))] <- NA
      if (any(is.na(lvls))) levels(tbl[[i]]) <- lvls
      
      lvls <- lvls[!is.na(lvls)]
      if (length(f <- which(lvls != trimws(lvls))) > 0) {
        field <- names(tbl)[[i]]
        cli_warn(
          c('i' = "Whitespace trimmed from {.field {field}} column {qty(length(f))} value{?s} {.val {lvls[f]}}.") )
        levels(tbl[[i]]) %<>% trimws()
      }
    }
  }
  
  attributes(tbl) <- attributes(tbl)[order(names(attributes(tbl)))]
  
  return (tbl)
}



require_package <- function (pkg, reason = 'for this command') {
  
  if (!nzchar(system.file(package = pkg)))
    package_missing(pkg = pkg, reason = reason)
  
  return (invisible(NULL))
}

package_missing <- function (pkg,  reason = 'for this command') {
  
  cli_abort(c(
    'x' = "The {.pkg {pkg}} R package is required {reason}.",
    'i' = "To install {.pkg {pkg}}, run:",
    '>' = if (!nzchar(system.file(package = 'pak'))) " {.run install.packages('pak')}",
    '>' = " {.run pak::pkg_install('{pkg}')}" ))
}

