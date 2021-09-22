# Helper functions called by the plot_*.r scripts.


#--------------------------------------------------------------
# Assign colors to categorical values.
#--------------------------------------------------------------
assign_colors <- function (vals, keys) {
  keys <- levels(keys)
  n    <- length(keys)
  
  if (is.null(vals)) {
    vals <- if (n <= 2) { # Colorblind friendly palette
      c('#00B9EB', '#ED5F16')
    } else if (n <= 8) {  # Colorblind friendly palette of 8 (jfly.iam.u-tokyo.ac.jp/color/)
      c("#0072B2", "#D55E00", "#CC79A7", "#009E73", "#F0E442", "#56B4E9", "#E69F00", "#999999")
    } else if (n <= 11) {   # Andrea's set of 11
      c('#66CDAA', '#FF8C69', '#8DB6CD', '#008B8B', '#FF6EB4', '#A2CD5A', '#FF6347', '#FFC125', 
        '#EEC591', '#BEBEBE', '#CD96CD')
    } else if (n <= 20) {   # Andrea's set of 20
      c('#8B4500', '#CD853F', '#FF1493', '#FFB5C5', '#CDC5BF', '#6C7B8B', '#B22222', '#FF0000', 
        '#FF7F24', '#FFD700', '#00FF00', '#2E8B57', '#00FFFF', '#63B8FF', '#0000FF', '#191970', 
        '#8B008B', '#A020F0', '#DA70D6', '#F08080')
    } else {                # Dan's set of 24
      c('#1C86EE', '#E31A1C', '#008B00', '#6A3D9A', '#A52A2A', '#FF7F00', '#FFD700', '#7EC0EE', 
        '#FB9A99', '#90EE90', '#CAB2D6', '#FDBF6F', '#B3B3B3', '#EEE685', '#B03060', '#FF83FA', 
        '#FF1493', '#0000FF', '#36648B', '#00CED1', '#00FF00', '#8B8B00', '#CDCD00', '#8B4500')
    }
  }
  
  assign_cleanup("colors", vals, keys)
}


#--------------------------------------------------------------
# Assign patterns to categorical values.
#--------------------------------------------------------------
assign_patterns <- function (vals, keys) {
  keys  <- levels(keys)
  
  dvals <- unique(c(
    'bricks', 'fishscales', 'right45', 'horizonal_saw', 
    'hs_cross', 'crosshatch45', 'crosshatch',
    ggpattern::magick_pattern_names ))
  
  if (is.null(vals)) {
    vals <- dvals
    
  } else if (is.null(names(vals)) && all(vals %in% keys)) {
    keys <- setNames(dvals[seq_along(vals)], vals)
  }
  
  assign_cleanup("patterns", vals, keys)
}


#--------------------------------------------------------------
# Assign shapes to categorical values.
#--------------------------------------------------------------
assign_shapes <- function (vals, keys) {
  keys <- levels(keys)
  n    <- length(keys)
  
  vals <- c(16, 17, 15, 3, 7, 8)
  if (n > 6)  vals <- c(0:14)
  if (n > 15) vals <- c(65:90, 97:122)
  if (n > 52) vals <- rep(vals, ceiling(n / 52))
  
  assign_cleanup("shapes", vals, keys)
}


#--------------------------------------------------------------
# Ensure correct length/keys for colors, patterns, & shapes
#--------------------------------------------------------------
assign_cleanup <- function (mode, vals, keys) {
  n <- length(keys)
  
  if (is.null(names(vals))) {
    if (length(vals) < n) vals <- rep_len(vals, n)
    if (length(vals) > n) vals <- vals[seq_len(n)]
    vals <- setNames(vals, keys)
    
  } else {
    missing <- setdiff(keys, names(vals))
    if (length(missing) > 3) missing <- c(head(missing, 3), "...")
    if (length(missing) > 0)
      stop("Missing ", mode, " assignments for: ", paste(sep=", ", missing))
    vals <- vals[keys]
  }
  
  return (vals)
}


#--------------------------------------------------------------
# Expand layer string / char vector to possible values
# x = "dr" or c("d", "bar") or more examples below
# choices = c(d = "dot", r = "bar")
# default = "dr"
#--------------------------------------------------------------
layer_match <- function (x, choices, default) {
  
  if (is.null(default))        default        <- choices[[1]]
  if (is.null(x))              x              <- default 
  # if (is.null(names(choices))) names(choices) <- unlist(substr(choices, 1, 1))
  
  x <- x[nchar(x) > 0]
  if (length(x) == 0) return (default)
  
  x    <- tolower(x)
  vals <- unname(choices)
  
  if (length(x) > 1) {
    x <- vals[pmatch(x, tolower(vals))] # c("d", "bar") => c("dot", "bar")
    
  } else if (nchar(x) == 1) {
    if (is.null(names(choices))) {
      x <- vals[pmatch(x, tolower(vals))]
    } else {
      x <- unname(choices[x])  # "r" => "bar"
    }
    
  } else if (is.na(pmatch(x, tolower(vals)))) {
    x <- strsplit(x, '')[[1]]
    
    if (is.null(names(choices))) {
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


#--------------------------------------------------------------
# Create the plot and add each layer with its arguments.
# Also attach a human-readable version of the plot command.
#--------------------------------------------------------------
layers_toPlot <- function (layers, dots) {
  
  stopifnot('ggplot' %in% names(layers))
  
  #--------------------------------------------------------------
  # 'layers' can draw on the layers listed here
  #--------------------------------------------------------------
  layerspecs <- list(
    'ggplot'     = list('.fn' = "ggplot",                    '.regex' = "^g(|gplot)\\."),
    'point'      = list('.fn' = "geom_point",                '.regex' = "^p(|oint)\\."),
    'violin'     = list('.fn' = "geom_violin",               '.regex' = "^v(|iolin)\\."),
    'box'        = list('.fn' = "geom_boxplot",              '.regex' = "^b(|ox)\\."),
    'crossbar'   = list('.fn' = "geom_crossbar",             '.regex' = "^c(|rossbar)\\."),
    'bar'        = list('.fn' = "geom_bar",                  '.regex' = "^(|ba)r\\."),
    'errorbar'   = list('.fn' = "geom_errorbar",             '.regex' = "^e(|rrorbar)\\."),
    'linerange'  = list('.fn' = "geom_linerange",            '.regex' = "^l(|inerange)\\."),
    'pointrange' = list('.fn' = "geom_pointrange",           '.regex' = "^p(|ointrange)\\."),
    'dot'        = list('.fn' = "geom_beeswarm",             '.regex' = "^d(|ot)\\."),
    'strip'      = list('.fn' = "geom_quasirandom",          '.regex' = "^s(|trip)\\."),
    'regression' = list('.fn' = "stat_smooth",               '.regex' = "^r(|egression)\\."),
    'stats_text' = list('.fn' = "geom_text",                 '.regex' = "^a(|nno)\\."),
    'brackets'   = list('.fn' = "geom_segment",              '.regex' = "^h(|line)\\."),
    'fill'       = list('.fn' = "scale_fill_manual",         '.regex' = "^fill\\."),
    'color'      = list('.fn' = "scale_color_manual",        '.regex' = "^color\\."),
    'shape'      = list('.fn' = "scale_shape_manual",        '.regex' = "^shape\\."),
    'pattern'    = list('.fn' = "scale_pattern_type_manual", '.regex' = "^pattern\\."),
    'facet'      = list('.fn' = "facet_wrap",                '.regex' = "^f(|acet)\\."),
    'labs'       = list('.fn' = "labs",                      '.regex' = "^labs\\."),
    'x_disc'     = list('.fn' = "scale_x_discrete",          '.regex' = "^x(|axis)\\."),
    'x_cont'     = list('.fn' = "scale_x_continuous",        '.regex' = "^x(|axis)\\."),
    'y_cont'     = list('.fn' = "scale_y_continuous",        '.regex' = "^y(|axis)\\."),
    'theme_bw'   = list('.fn' = "theme_bw",                  '.regex' = "^(theme_|)bw\\."),
    'theme'      = list('.fn' = "theme",                     '.regex' = "^t(|heme)\\.") )
  
  
  #--------------------------------------------------------------
  # Other common settings
  #--------------------------------------------------------------
  layerspecs[['stats_text']][['show.legend']] <- FALSE
  layerspecs[['errorbar']][['width']] <- 0.5
  layerspecs[['crossbar']][['width']] <- 0.5
  
  layerspecs[['stats_text']] %<>% c(list(
    'color'   = "black",
    'mapping' = list(x=".x", y=".y", label=".label") ))
  layerspecs[['brackets']] %<>% c(list(
    'color'   = "black",
    'mapping' = list(x=".x", y=".y", xend=".xend", yend=".yend") ))
  
  for (i in c('errorbar', 'crossbar', 'linerange', 'pointrange')) {
    layerspecs[[i]][['mapping']] <- list(y=".y", ymin=".ymin", ymax=".ymax")
  }
  
  
  #--------------------------------------------------------------
  # Only subset the data when multiple data layers are present
  #--------------------------------------------------------------
  if (isTRUE(length(unique(layers[['ggplot']][['data']][['.src']])) > 1)) {
    layerspecs[['stats_text']][['data']] <- ~ subset(., .src == "stats")
    layerspecs[['brackets']][['data']]   <- ~ subset(., .src == "brackets")
    
    for (i in c('point', 'bar', 'violin', 'strip', 'box', 'dot', 'regression')) {
      layerspecs[[i]][['data']] <- ~ subset(., .src == "points")
    }
    
    for (i in c('errorbar', 'crossbar', 'linerange', 'pointrange')) {
      layerspecs[[i]][['data']] <- ~ subset(., .src == "vline")
    }
  }
  
  
  #--------------------------------------------------------------
  # Special cases for ggbeeswarm and ggpattern functions
  #--------------------------------------------------------------
  patterned <- isTRUE("pattern" %in% names(layers))
  for (layer in names(layers)) {
    
    if (!layer %in% names(layerspecs)) {
      stopifnot(isTRUE(nzchar(layers[[layer]][['.fn']])))
      stopifnot(is.function(layers[[layer]][['.fun']]))
    }
    
    if (patterned && layer %in% c("bar", "box", "crossbar", "violin")) {
      layerspecs[[layer]][['.fn']] %<>% paste0(., "_pattern")
    }
    
    if (patterned && layer %in% c("fill", "color", "shape")) {
      layerspecs[[layer]][['.fn']] <- paste0("scale_pattern_", layer, "_manual")
    }
    
    pkg <- "ggplot2"
    if (grepl("_beeswarm",    layerspecs[[layer]][['.fn']], fixed = TRUE)) pkg <- "ggbeeswarm"
    if (grepl("_quasirandom", layerspecs[[layer]][['.fn']], fixed = TRUE)) pkg <- "ggbeeswarm"
    if (grepl("_pattern",     layerspecs[[layer]][['.fn']], fixed = TRUE)) pkg <- "ggpattern"
    
    layerspecs[[layer]][['.fun']] <- do.call(`::`, list(pkg, layerspecs[[layer]][['.fn']]))
    if (pkg != "ggplot2") layerspecs[[layer]][['.fn']] %<>% paste0(pkg, "::", .)
    
    # Merge custom and default top-level parameters
    for (i in setdiff(names(layerspecs[[layer]]), names(layers[[layer]])))
      layers[[layer]][[i]] <- layerspecs[[layer]][[i]]
    
    # Merge custom and default mapping parameters
    for (i in setdiff(names(layerspecs[[layer]][['mapping']]), names(layers[[layer]][['mapping']])))
      layers[[layer]][['mapping']][[i]] <- layerspecs[[layer]][['mapping']][[i]]
  }
  
  
  
  #--------------------------------------------------------------
  # Standardize the list order: ggplot() first, theme() last, etc
  #--------------------------------------------------------------
  layers <- layers[c(
    intersect(names(layerspecs), names(layers)),
    setdiff(names(layers), names(layerspecs))
  )]
  
  
  p    <- NULL
  cmds <- c("library(ggplot2)")
  
  for (layer in names(layers)) {
    
    args <- layers[[layer]]
    fn   <- args[['.fn']]
    fun  <- args[['.fun']]
    regx <- args[['.regex']]
    args <- args[grep("^\\.", names(args), invert = TRUE)]
    
    
    # Unprefixed dot arguments, e.g. 'scales'="free_x"
    #--------------------------------------------------------------
    for (i in intersect(names(dots), formalArgs(fun)))
      args[[i]] <- dots[[i]]
    
    
    # Prefixed dot arguments, e.g. 'facet.scales'="free_x"
    #--------------------------------------------------------------
    for (i in grep(regx, names(dots), value = TRUE))
      args[[sub(regx, "", i, perl = TRUE)]] <- dots[[i]]
    
    
    # Don't specify an argument if it's already the default
    #--------------------------------------------------------------
    defaults <- formals(fun)
    for (i in intersect(names(defaults), names(args)))
      if (identical(defaults[[i]], args[[i]]))
        args[[i]] <- NULL
    
    
    # Create the aes object for mapping=
    #--------------------------------------------------------------
    if ('mapping' %in% names(args))
      args[['mapping']] <- do.call(
        what  = aes_string, 
        args  = args[['mapping']], 
        quote = TRUE )
    
    
    # Rewrite layer with the updated args
    #--------------------------------------------------------------
    layers[[layer]] <- args
    
    
    # Skip theme() unless it has arguments
    #--------------------------------------------------------------
    if (fn == "theme" && length(args) == 0) next
    
    
    # Show ggplot() layer as "ggplot(data)", rest more verbosely
    #--------------------------------------------------------------
    if (fn == "ggplot") {
      p <- do.call(fun, args) 
      cmds <- sprintf("%s(%s)", fn, as.args(args, fun=fun))
    } else {
      p    <- p + do.call(fun, args)
      cmds <- c(cmds, sprintf("%s(%s)", fn, as.args(args, indent = 4, fun=fun)))
    }
    
  }
  
  #--------------------------------------------------------------
  # Attach the number of facet rows and cols as plot attributes
  #--------------------------------------------------------------
  if ('nfacets' %in% names(attributes(dots))) {
    rc <- ggplot2::wrap_dims(
      n    = attr(dots, 'nfacets'), 
      nrow = layers[['facet']][['nrow']], 
      ncol = layers[['facet']][['ncol']] )
    attr(dots, 'facet.nrow') <- min(rc[[1]], attr(dots, 'nfacets'))
    attr(dots, 'facet.ncol') <- min(rc[[2]], attr(dots, 'nfacets'))
  }
  attr(p, 'facet.nrow') <- attr(dots, 'facet.nrow') %||% 1
  attr(p, 'facet.ncol') <- attr(dots, 'facet.ncol') %||% 1
  
  
  attr(p, 'cmd')  <- paste(collapse=" +\n  ", cmds)
  attr(p, 'data') <- layers[['ggplot']][['data']]
  
  return (p)
}


#------------------------------------------------------------------
# Convert a list of arguments to character strings.
# When indent > 0, produces a multi-line string.
# When fun is a function, puts the arguments in the expected order.
#------------------------------------------------------------------
as.args <- function (args = list(), indent = 0, fun = NULL) {
  
  stopifnot(is.list(args))
  stopifnot(is.numeric(indent))
  
  # Right-pad parameter names.
  fmt <- "%s = %s"
  if (isTRUE(indent > 0 && length(args) > 1))
    fmt <- paste(collapse="", rep(" ", indent)) %>%
    paste0("\n", ., "%-", max(nchar(names(args))), "s = %s")
  
  # Re-arrange parameter order; omit names where implicitly known.
  if (is.function(fun)) {
    f_args <- formalArgs(fun)
    args   <- args[c(intersect(f_args, names(args)), sort(setdiff(names(args), f_args)))]
    if (isTRUE(indent == 0))
      for (i in seq_along(args))
        if (names(args)[[i]] == f_args[[i]]) names(args)[i] <- "" else break
  }
  
  # Convert arguments to `eval`-able string representations
  strs <- c()
  for (i in seq_along(args)) {
    
    key <- names(args)[[i]]
    val <- args[[i]]
    
    val <- if (is.null(val))                    { "NULL"
    } else if (!is.null(attr(val, 'display')))  { attr(val, 'display') 
    } else if (is.character(val))               { glue::double_quote(val) 
    } else if (is.logical(val))                 { as.character(val)
    } else if (is.numeric(val))                 { as.character(val)
    } else if (is(val, 'BIOM'))                 { "biom"
    } else if (is.data.frame(val))              { "data"
    } else if (is(val, 'formula'))              { capture.output(val)[[1]]
    } else if (is.function(val))                { fun_toString(val)
    } else if (is(val, 'uneval'))               { aes_toString(val)
    } else if (is.factor(val))                  { as.character(val)
    } else                                      { capture.output(val) }
    
    
    if (isTRUE(is.character(val) && !is.null(names(val))))
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


#------------------------------------------------------------------
# Find a function's name
#------------------------------------------------------------------
fun_toString <- function (x) {
  
  pkg <- environment(x)[['.packageName']]
  if (!is.null(pkg))
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


#------------------------------------------------------------------
# Convert an aes object to a string
#------------------------------------------------------------------
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
    val <- if (is(val, 'formula')) { capture.output(as.name(all.vars(val)))
    } else                         { glue::double_quote(val) }
    
    key %<>% sub(pattern = "colour", replacement = "color")
    key <- capture.output(as.name(key))
    
    results %<>% c(sprintf("%s = %s", key, val))
  }
  
  return (sprintf("aes(%s)", paste(collapse = ", ", results)))
}


#------------------------------------------------------------------
# rbind(), but add/rearrange columns as needed
#------------------------------------------------------------------
append_df <- function (x, y) {
  xy <- unique(c(names(x), names(y)))
  for (i in setdiff(xy, names(x))) x[[i]] <- NA
  for (i in setdiff(xy, names(y))) y[[i]] <- NA
  rbind(x[,xy,drop=F], y[,xy,drop=F])
}


#------------------------------------------------------------------
# Explicitly define the code to be displayed in cmd
#------------------------------------------------------------------
as.cmd <- function (expr, env=NULL) {
  if (is.null(env)) {
    cmd <- capture.output(substitute(expr))
  } else {
    cmd <- do.call(substitute, list(expr=substitute(expr), env=env))
    cmd <- capture.output(cmd)
  }
  
  if (length(cmd) > 1)
    cmd <- paste(trimws(cmd), collapse = " ")
  
  structure(expr, 'display' = cmd)
}


#------------------------------------------------------------------
# Identify and remove rows of df with bad values
#------------------------------------------------------------------
finite_check <- function (df, col=".y", metric=NULL) {
  
  if (all(is.finite(df[[col]])))
    return (NULL)
  
  bad <- which(!is.finite(df[[col]]))
  n   <- length(bad)
  
  if (is.null(metric))
    if (".metric" %in% names(df))
      if (length(unique(df[['.metric']])) == 1)
        metric <- df[1,'.metric']
  
  metric <- ifelse(is.null(metric), "", paste0(metric, " "))
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

