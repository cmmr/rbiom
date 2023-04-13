# Helper functions called by the plot_*.r scripts.

#____________________________________________________________________
# Assign RHS to LHS if LHS is NULL.
#____________________________________________________________________
`%||=%` <- function(.lhs, .rhs) {
  .lhs_quo <- rlang::enquo(.lhs)
  .lhs_str <- rlang::quo_name(.lhs_quo)
  .lhs_env <- rlang::get_env(.lhs_quo)
  
  if (!exists(.lhs_str, envir = rlang::get_env(.lhs_quo)))
    .lhs <- NULL
  
  assign(
    x     = .lhs_str, 
    value = if (is_null(.lhs)) .rhs else .lhs, 
    envir = .lhs_env )
}


#____________________________________________________________________
# Expand layer string / char vector to possible values
# x = "dr" or c("d", "bar") or more examples below
# choices = c(d = "dot", r = "bar")
# default = "dr"
#____________________________________________________________________
layer_match <- function (x, choices, default) {
  
  if (is_null(default))        default        <- choices[[1]]
  if (is_null(x))              x              <- default 
  # if (is_null(names(choices))) names(choices) <- unlist(substr(choices, 1, 1))
  
  x <- x[nchar(x) > 0]
  if (length(x) == 0) return (default)
  if (all(nchar(x) == 1)) x <- paste(collapse = "", x)
  
  x    <- tolower(x)
  vals <- unname(choices)
  
  if (length(x) > 1) {
    x <- vals[pmatch(x, tolower(vals))] # c("d", "bar") => c("dot", "bar")
    
  } else if (nchar(x) == 1) {
    if (is_null(names(choices))) {
      x <- vals[pmatch(x, tolower(vals))]
    } else {
      x <- unname(choices[x])  # "r" => "bar"
    }
    
  } else if (is.na(pmatch(x, tolower(vals)))) {
    x <- strsplit(x, '')[[1]]
    
    if (is_null(names(choices))) {
      x <- vals[pmatch(x, tolower(vals))]
    } else {
      x <- unname(choices[x]) # "dr" => c("dot", "bar")
    }
    
  } else {
    x <- vals[pmatch(x, tolower(vals))] # "do" => "dot"
  }
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return (default)
  
  return (x)
}


#____________________________________________________________________
# Convert a list of arguments to character strings.
# When indent > 0, produces a multi-line string.
# When fun is a function, puts the arguments in the expected order.
#____________________________________________________________________
as.args <- function (args = list(), indent = 0, fun = NULL) {
  
  stopifnot(is_list(args))
  stopifnot(is_integerish(indent))
  
  
  # Discard arguments with `display = FALSE` attribute.
  for (i in rev(seq_along(args)))
    if (isFALSE(attr(args[[i]], 'display', exact = TRUE)))
      args[[i]] <- NULL
  
  
  # Re-arrange parameter order; omit names where implicitly known; ignore if default.
  if (is.function(fun) && !is_null(names(args))) {
    f_args <- formals(fun)
    
    for (i in names(args))
      if (hasName(f_args, i) && identical(args[[i]], f_args[[i]]))
        args[[i]] <- NULL
    
    f_args <- names(f_args)
    args   <- args[c(intersect(f_args, names(args)), sort(setdiff(names(args), f_args)))]
    if (isTRUE(indent == 0))
      for (i in seq_along(args))
        if (names(args)[[i]] == f_args[[i]]) names(args)[i] <- "" else break
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
    
    val <- if (is_null(val))        { "NULL"
    } else if (!is_null(display))   { display
    } else if (is_character(val))   { glue::double_quote(val) 
    } else if (is_logical(val))     { as.character(val) %>% setNames(names(val))
    } else if (is.numeric(val))     { as.character(val) %>% setNames(names(val))
    } else if (is(val, 'quosures')) { as.character(val)
    } else if (is(val, 'BIOM'))     { "biom"
    } else if (is.data.frame(val))  { "data"
    } else if (is(val, 'formula'))  { capture.output(val)[[1]]
    } else if (is.function(val))    { fun_toString(val)
    } else if (is(val, 'uneval'))   { aes_toString(val)
    } else if (is.factor(val))      { as.character(val) %>% setNames(names(val))
    } else if (is_list(val))        { paste0("list(", as.args(val), ")")
    } else                          { capture.output(val) }
    
    
    if (isTRUE(is_character(val) && !is_null(names(val))))
      val <- paste(glue::single_quote(names(val)), "=", unname(val))
    
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
      if (identical(x, do.call(`::`, list(pkg, fn)))) {
        if (!pkg %in% getOption("defaultPackages"))
          fn <- sprintf("%s::%s", pkg, fn)
        return (fn)
      }
  
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
    val <- if (is(val, 'quosure')) { rlang::as_label(val)
    } else if (is(val, 'formula')) { capture.output(as.name(all.vars(val)))
    } else if (is.logical(val))    { as.character(val)
    } else if (is.numeric(val))    { as.character(val)
    } else                         { glue::double_quote(val) }
    
    key %<>% sub(pattern = "colour", replacement = "color")
    key <- capture.output(as.name(key))
    
    results %<>% c(sprintf("%s = %s", key, val))
  }
  
  return (sprintf("aes(%s)", paste(collapse = ", ", results)))
}


#____________________________________________________________________
# rbind(), but add/rearrange columns as needed
#____________________________________________________________________
append_df <- function (...) {
  
  Reduce(
    x = Filter(f = Negate(is.null), x = list(...)),
    f = function (x, y) {
      xy <- unique(c(names(x), names(y)))
      for (i in setdiff(xy, names(x))) x[[i]] <- NA
      for (i in setdiff(xy, names(y))) y[[i]] <- NA
      rbind(x[,xy,drop=F], y[,xy,drop=F])
    })
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
  cols <- unlist(list(...))
  if (is_null(cols) || is_null(df)) return (NULL)
  df[,intersect(colnames(df), cols),drop=FALSE]
}
rename_cols <- function(df, ...) {
  vals <- list(...)
  for (i in names(vals))
    names(df)[which(names(df) == i)] <- vals[[i]]
  return (df)
}


#____________________________________________________________________
# Rename common data frame column names, avoiding conflicts
#____________________________________________________________________
soft_rename <- function (df, safe=FALSE) {
  
  if (safe)
    return (df)
  
  map <- c(
    ".test"   = "Test",        ".metric" = "Metric",
    ".sample" = "Sample",      ".taxa"   = "Taxa", 
    ".p.val"  = "P-Value",     ".adj.p"  = "Adjusted P",
    ".f.stat" = "F-Statistic", ".r.sqr"  = "R-Squared",
    ".x"      = "x",           ".ori.1"  = "ori.1", 
    ".y"      = "y",           ".ori.2"  = "ori.2", 
    ".axis.1" = "Axis 1", 
    ".axis.2" = "Axis 2",
    ".value"  = "value" )
  
  for (i in names(map))
    if (i %in% names(df) && !map[[i]] %in% names(df))
      names(df)[which(names(df) == i)] <- map[[i]]
  
  return (df)
}


#____________________________________________________________________
# Decorate column names for plyr
#____________________________________________________________________
ply_cols <- function (cols) {
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
      glue::glue_collapse(
        x = paste(df[bad,'.sample'], "=", df[bad,col]), 
        width = 100, sep = ", ", last = ", and " ))
  } else {
    msg %<>% paste0(
      glue::glue_collapse(
        x = as.character(unique(df[bad,col])), 
        width = 100, sep = ", ", last = ", and " ))
  }
  
  list('bad' = bad, 'msg' = msg)
}





#____________________________________________________________________
# Log scale labels
#____________________________________________________________________
si_units <- as.cmd(scales::label_number(scale_cut = scales::cut_si("")))

loglabels <- function (values) {
  
  hi <- ceiling(log10(max(force(values))))
  list(
    'breaks'       = as.cmd(10 ** (0:hi),                    env = list(hi = hi)),
    'minor_breaks' = as.cmd(as.vector(10 ** (0:hi) %o% 2:9), env = list(hi = hi - 1)),
    'labels'       = si_units )
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
# Test if ALL arguments are NULL
#____________________________________________________________________
all.null <- function (...) {
  for (i in list(...))
    if (!is_null(i))
      return (FALSE)
  return (TRUE)
}

#____________________________________________________________________
# Test if ANY arguments are NULL
#____________________________________________________________________
any.null <- function (...) {
  for (i in list(...))
    if (is_null(i))
      return (TRUE)
  return (FALSE)
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


