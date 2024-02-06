
# Helper functions called by the plot_*.r scripts.



#____________________________________________________________________
# Stop with an error message if variable isn't an rbiom object.
#____________________________________________________________________
must_be_rbiom <- function (biom) {
  if (!(is(biom, 'rbiom') && is(biom, 'R6')))
    cli_abort(('x' = "`biom` must be an rbiom object, not {.type {biom}}."))
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
# Expand layer string / char vector to possible values
# x = "dr" or c("d", "bar") or more examples below
# choices = c(d = "dot", r = "bar")
# default = "dr"
#____________________________________________________________________
layer_match <- function (x, choices, default) {
  
  result <- c()
  
  for (i in tolower(x[nchar(x) > 0]))
    result <- c(result, local({
      
      if (i %in% names(choices)) return (choices[[i]])
      
      pm <- pmatch(i, choices)
      ii <- strsplit(i, '')[[1]]
      
      if (length(x) > 1 && !is.na(pm)) return (choices[[pm]])
      if (all(ii %in% names(choices))) return (choices[ii])
      if (!is.na(pm))                  return (choices[[pm]])
      return (choices[intersect(ii, names(choices))])
    }))
  
  result <- result[!is.na(result)]
  
  if (length(result) == 0)
    result <- choices[strsplit(default, '')[[1]]]
  
  return (result)
}




#____________________________________________________________________
# Convert to tibble and add 'rbiom_tbl' class.
#____________________________________________________________________
as_rbiom_tbl <- function (df) {
  add_class(as_tibble(df), 'rbiom_tbl')
}



#____________________________________________________________________
# Create an environment of evaluated variables in caller's env.
#____________________________________________________________________
eval_envir <- function (env, ...) {
  
  dots <- list(...)
  
  
  #____________________________________________________________________
  # Evaluate all arguments into a clean env.
  #____________________________________________________________________
  params <- lapply(as.list(env), eval) %>%
    rlang::new_environment(parent = rlang::ns_env('rbiom'))
  
  
  #____________________________________________________________________
  # Move some parameters into dots.
  #____________________________________________________________________
  if (hasName(params, 'y.trans')) {
    dots[['y.trans']] <- params[['y.trans']]
    rlang::env_unbind(params, 'y.trans')
  }
  
  
  #____________________________________________________________________
  # Check if obviously wrong parameters are going into dots.
  #____________________________________________________________________
  if (length(dots) > 0) {
    
    invalid <- paste(collapse = ", ", c(
      names(dots)[startsWith(names(dots), '.')],
      intersect(
        x = names(dots),
        y = c(
          'x', 'color.by', 'shape.by', 'facet.by', 
          'pattern.by', 'label.by', 'order.by', 'stat.by', 'limit.by',
          'adiv', 'bdiv', 'taxa', 'rank', 'weighted', 
          'within', 'between', 'tree', 'params' ))))
    
    if (nzchar(invalid)) {
      fn <- deparse(rlang::caller_call()[[1]])
      cli_abort(c(x = "{fn} does not accept parameter(s): {invalid}"))
    }
    
  }
  params$.dots <- lapply(dots, eval)
  
  
  return (params)
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
# Run a command and add attr(,cmd) and attr(,display) to result.
#____________________________________________________________________
run.cmd <- function (f, args, hist=NULL, lhs=NULL, display=lhs, envir = parent.frame()) {
  
  fn  <- as.character(substitute(f))
  fun <- f
  if (is_scalar_character(f)) {
    fn  <- f
    x   <- strsplit(f, "::", fixed = TRUE)[[1]]
    fun <- if (length(x) == 1) get(x[[1]]) else getFromNamespace(x[[2]], ns = x[[1]])
  }
  
  cmd <- sprintf("%s(%s)", fn, as.args(args, fun=fun))
  if (is.null(display)) display <- if.null(lhs, cmd)
  if (!is.null(lhs))    cmd <- sprintf("%s <- %s", lhs, cmd)
  if (!is.null(hist))   cmd <- paste(collapse = "\n", c(attr(hist, 'code'), cmd))
  
  return (aa(do.call(fun, args), display = display, code = cmd))
}


#____________________________________________________________________
# Take attr(,'code') from `hist` and append `lhs <- fn(args)`.
#____________________________________________________________________
append.cmd <- function (x, fn, args, lhs = NULL, hist = NULL, display = lhs, fun = get(fn)) {
  
  cmd <- sprintf("%s(%s)", fn, as.args(args, fun=fun))
  
  if (is.null(display)) display <- if.null(lhs, cmd)
  if (!is.null(lhs))    cmd <- sprintf("%s <- %s", lhs, cmd)
  if (!is.null(hist))   cmd <- paste(collapse = "\n", c(attr(hist, 'code'), cmd))
  
  return (aa(x, display = display, code = cmd))
}



#____________________________________________________________________
# Convert a list of arguments to character strings.
# When indent > 0, produces a multi-line string.
# When fun is a function, puts the arguments in the expected order.
#____________________________________________________________________
as.args <- function (args = list(), indent = 0, fun = NULL) {
  
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
    if (isTRUE(indent == 0))
      for (i in seq_along(args))
        if (eq(names(args)[[i]], f_args[[i]])) names(args)[i] <- "" else break
  }
  
  # Right-pad parameter names.
  fmt <- "%s = %s"
  if (isTRUE(indent > 0 && length(args) > 1 && !is_null(names(args))))
    fmt <- paste(collapse="", rep(" ", indent)) %>%
    paste0("\n", ., "%-", max(nchar(names(args))), "s = %s")
  
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
    } else if (is(val, 'quosures'))   { as.character(val)
    } else if (is(val, 'rbiom'))      { "biom"
    } else if (is.data.frame(val))    { "data"
    } else if (is(val, 'formula'))    { format(val)
    } else if (is.function(val))      { fun_toString(val)
    } else if (is(val, 'uneval'))     { aes_toString(val)
    } else if (is.factor(val))        { as.character(val) %>% setNames(names(val))
    } else if (is_list(val))          { paste0("list(", as.args(val), ")")
    } else                            { capture.output(val) }
    
    
    if (isTRUE(is_character(val) && !is_null(names(val))))
      val <- paste(single_quote(names(val)), "=", unname(val))
    
    if (length(val) > 1)
      val <- paste0("c(", paste(collapse=", ", val), ")")
    
    if (nzchar(key)) {
      key <- capture.output(as.name(key))
      val <- sprintf(fmt = fmt, key, val)
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
    "colour", "fill", "shape", "group", 
    "pattern_type", "pattern_fill", "label" )
  keys <- c(intersect(sort, keys), setdiff(keys, sort))
  
  results <- c()
  for (key in keys) {
    
    val <- x[[key]]
    val <- if (is(val, 'quosure')) { capture.output(rlang::quo_get_expr(val))
    } else if (is(val, 'formula')) { capture.output(as.name(all.vars(val)))
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
# Example: fmt("dm <- %s", bdiv_distmat, biom, weighted, tree)
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
# rbind(), but add/rearrange columns as needed
#____________________________________________________________________
append_df <- function (...) {
  
  Reduce(
    x = Filter(f = Negate(is.null), x = list(...)),
    f = function (x, y) ({
      xy <- unique(c(names(x), names(y)))
      for (i in setdiff(xy, names(x))) x[[i]] <- NA
      for (i in setdiff(xy, names(y))) y[[i]] <- NA
      rbind(x[,xy,drop=F], y[,xy,drop=F])
    }))
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
drop_empty <- function (df) {
  if (is_null(df))   return (NULL)
  if (nrow(df) == 0) return (df)
  df[,apply(df, 2L, function (x) !all(is.na(x))),drop=FALSE]
}
keep_cols <- function (df, ...) {
  cols <- unique(unlist(list(...)))
  if (is_null(cols) || is_null(df)) return (NULL)
  df[,intersect(cols, colnames(df)),drop=FALSE]
}
rename_cols <- function(df, ...) {
  
  vals <- list(...)
  if (length(vals) == 1) vals <- as.list(vals[[1]])
  
  x <- attr(df, 'names', exact = TRUE)
  for (i in names(vals))
    if (i %in% x)
      x[which(x == i)] <- vals[[i]]
  
  attr(df, 'names') <- x
  return (df)
}
rename_response <- function(df, new) {
  old <- attr(df, 'response', exact = TRUE)
  stopifnot(is_scalar_character(old) && !is_na(old))
  stopifnot(is_string(old, colnames(old)))
  names(df)[which(names(df) == old)] <- new
  attr(df, 'response') <- new
  return (df)
}


#____________________________________________________________________
# Decorate column names for plyr
#____________________________________________________________________
ply_cols <- function (cols) {
  
  if (is(cols, 'quoted')) return (cols)
  
  vars <- NULL
  for (col in cols)
    vars <- c(plyr::as.quoted(as.name(col)), vars)
  
  return (vars)
}


#____________________________________________________________________
# Insert a new named column after a specified column index/name
#____________________________________________________________________
insert_col <- function (df, after=0, col, val) {
  
  if (is_null(df))                  return (NULL)
  if (is_null(col) || is_null(val)) return (df)
  
  df <- cbind(df, val)
  colnames(df)[ncol(df)] <- col
  
  
  after <- after %||% 0
  if (is.character(after))
    after <- which(colnames(df) == after) %>% head(1) %||% 0
  if (after > ncol(df) - 1)
    after <- ncol(df)
  
  pos <- unique(c(seq_len(after), ncol(df), seq_len(ncol(df))))
  df  <- df[,pos,drop=F]
  
  return (df)
}


#____________________________________________________________________
# Identify and remove rows of df with bad values
#____________________________________________________________________
finite_check <- function (df, col=".y", metric=NULL) {
  
  if (all(is.finite(df[[col]])))
    return (NULL)
  
  bad <- which(!is.finite(df[[col]]))
  n   <- length(bad)
  
  if (is_null(metric))
    if (".metric" %in% names(df))
      if (length(unique(df[['.metric']])) == 1)
        metric <- df[1,'.metric']
  
  metric <- ifelse(is_null(metric), "", paste0(metric, " "))
  msg    <- ifelse(n == 1, 
                   paste0("One sample had a non-finite ", metric, "value and was excluded from this plot:\n"),
                   paste0(n, " samples had non-finite ", metric, "values and were excluded from this plot:\n") )
  
  
  if ('.sample' %in% names(df)) {
    msg %<>% paste0(
      glue_collapse(
        x = paste(df[bad,'.sample'], "=", df[bad,col]), 
        width = 100, sep = ", ", last = ", and " ))
  } else {
    msg %<>% paste0(
      glue_collapse(
        x = as.character(unique(df[bad,col])), 
        width = 100, sep = ", ", last = ", and " ))
  }
  
  list('bad' = bad, 'msg' = msg)
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
# Turn unquoted barewords into a character vector.
#____________________________________________________________________
qw <- function (...) {
  all.vars(match.call())
}

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
  capture.output(as.name(nm))
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
# String describing object's class.
#____________________________________________________________________
of_class <- function (x) {
  paste0(" of class ", paste(collapse = " / ", class(x)))
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
# Shorthand for replacing NAs with something else.
#____________________________________________________________________
if.na <- function (x, replacement) {
  if (!is.na(x))                 return (x)
  if (!is.function(replacement)) return (replacement)
  return (replacement())
}

#____________________________________________________________________
# Shorthand for replacing NULLs with something else.
#____________________________________________________________________
if.null <- function (x, replacement) {
  if (!is_null(x))               return (x)
  if (!is.function(replacement)) return (replacement)
  return (replacement())
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
  
  stopifnot(is.function(fun))
  
  if (is.environment(params))
    params <- as.list(params, all.names = TRUE)
  
  stopifnot(is.list(params))
  
  if (is.list(params$.dots))
    params <- c(params[names(params) != ".dots"], params$.dots)
  
  fetch <- intersect(formalArgs(fun), names(params))
  
  return (params[fetch])
}

