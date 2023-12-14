

validate_adiv <- function (var = "adiv", env = parent.frame(), ...) {
  choices <- c("OTUs", "Shannon", "Chao1", "Simpson", "InvSimpson")
  validate_var_choices(var, choices, env, all_option = ".all", ...)
}

validate_bdiv <- function (var = "bdiv", env = parent.frame(), ...) {
  choices <- c("UniFrac", "Jaccard", "Bray-Curtis", "Manhattan", "Euclidean")
  validate_var_choices(var, choices, env, all_option = ".all", ...)
}

validate_ord <- function (var = "ord", env = parent.frame(), ...) {
  choices <- c("PCoA", "tSNE", "NMDS", "UMAP")
  validate_var_choices(var, choices, env, all_option = ".all", ...)
}

validate_dist <- function (var = "dist", env = parent.frame(), ...) {
  choices <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  validate_var_choices(var, choices, env, all_option = ".all", ...)
}

validate_clust <- function (var = "clust", env = parent.frame(), ...) {
  choices <- c("average", "ward", "mcquitty", "single", "median", "complete", "centroid")
  validate_var_choices(var, choices, env, all_option = ".all", ...)
}

validate_unc <- function (var = "unc", env = parent.frame(), ...) {
  choices <- c("asis", "singly", "grouped", "drop")
  validate_var_choices(var, choices, env, ...)
}



validate_tree <- function (var = "tree", env = parent.frame(), evar = var, null_ok = FALSE) {
  
  x <- get(var, pos = env, inherits = FALSE)
  
  if (is.null(x) && null_ok) return (invisible(NULL))
  if (is(x, 'phylo'))        return (invisible(NULL))
  
  if (is_scalar_character(x)) {
    x <- read_tree(x)
    assign(var, x, pos = env)
    return (invisible(NULL))
  }
  
  stop("`", evar, "` must be a phylo-class object, such as from read_tree().")
}





validate_taxa <- function (var = "taxa", env = parent.frame(), evar = var, null_ok = FALSE) {
  
  x <- get(var, pos = env, inherits = FALSE)
  
  if (is.null(x) && null_ok) return (invisible(NULL))
  if (is.null(x))             stop ("`", evar, "` can't be NULL.")
  if (anyNA(x))               stop ("`", evar, "` can't be NA.")
  
  if (!(is.numeric(x) || is.character(x)))
    stop ("`", evar, "` must be one number or a character vector.")
  
  if (is.numeric(x) && length(x) > 1)
    stop ("`", evar, "` must be one number or a character vector.")
  
  if (is.numeric(x) && x <= 0)
    stop ("When numeric, `", evar, "` must be greater than zero.")
  
  return (invisible(NULL))
}





validate_rank <- function (var = "rank", env = parent.frame(), evar = var, null_ok = FALSE, ...) {
  
  #   0     1      2      3     4     5      6
  # .otu Kingdom Phylum Class Order Family Genus
  #   0    -6     -5     -4    -3    -2     -1
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  x       <- ifelse(tolower(x) == "otu", ".otu", x)
  biom    <- get("biom", pos = env, inherits = FALSE)
  choices <- biom$ranks
  
  
  if (is_integerish(x)) {
    x <- ifelse(x >= 0, choices[x + 1], rev(choices)[abs(x)])
    if (anyNA(x)) {
      n <- length(choices)
      stop ("When numeric, `", evar, "` must be between -", n, " and ", n, ".")
    }
  }
  
  validate_var_choices('x', choices, null_ok = null_ok, evar = var, ...)
  assign(var, x, pos = env)
  
  return (invisible(NULL))
}



validate_meta <- function (var = "meta", env = parent.frame(), evar = var, null_ok = FALSE, max = 1, cmp = FALSE, col_type = NULL) {
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  biom    <- get("biom", pos = env, inherits = FALSE)
  choices <- biom$fields
  is_num  <- sapply(USE.NAMES = TRUE, choices, function (i) is.numeric(biom$metadata[[i]]))
  
  
  if (is_scalar_logical(x)) {
    
    warning(evar, " should be a character vector, not ", x, ".")
    
    if (isTRUE(x))  x <- choices
    if (isFALSE(x)) x <- NULL
  
    # Exit early if TRUE/FALSE resolved to NULL
    if (is.null(x) && null_ok) {
      assign(var, x, pos = env)
      return (invisible(NULL))
    }
  }
  
  
  if (eq(x, '.all')) {
    x <- choices[-1]
    
    # Exit early if '.all' resolved to NULL
    if (is.null(x) && null_ok) {
      assign(var, x, pos = env)
      return (invisible(NULL))
    }
  }
  
  
  if (is_null(x))       stop ("`", evar, "` cannot be NULL.")
  if (!is_character(x)) stop ("Invalid argument to `", evar, "`", of_class(x))
  if (length(x) > max)  stop ("`", evar, "` cannot be length ", length(x), ".")
  if (anyNA(x))         stop ("`", evar, "` cannot be NA.")
  x <- trimws(x)
  if (!all(nzchar(x)))  stop ("`", evar, "` cannot be ''.")
  
  
  within  <- NULL
  between <- NULL
  revsort <- NULL
  
  for (i in seq_along(x)) {
    
    xi <- x[[i]]
    
    
    # The prefixes '==' and '!=' are used in bdiv ops to limit comparisons.
    if (startsWith(xi, "==") || startsWith(xi, "!=")) {
      col_name <- substr(xi, 3, nchar(xi))
      if (!isTRUE(cmp))       stop ("Can't use prefixes, e.g. '", xi, "', in this function.")
      if (is_num[[col_name]]) stop ("Invalid prefix for numeric metadata field: '", xi, "'.")
      if (startsWith(xi, "==")) within  %<>% c(i)
      if (startsWith(xi, "!=")) between %<>% c(i)
      xi <- col_name
    }
    
    
    # A prefix of '-' to indicates reverse sorting.
    if (startsWith(xi, "-") && !xi %in% choices) {
      col_name <- substr(xi, 2, nchar(xi))
      if (!is_num[[col_name]]) stop ("Invalid prefix for categorical metadata field: '", xi, "'.")
      revsort %<>% c(i)
      xi <- col_name
    }
    
    
    # Check for valid name and (optionally) type.
    if (!xi %in% choices)
      validate_var_choices('xi', choices, evar = evar)
    
    if (eq(col_type, "cat") && is_num[[xi]])
      stop ("`", evar, "` field '", xi, "' is not categorical.")
    
    if (eq(col_type, "num") && !is_num[[xi]])
      stop ("`", evar, "` field '", xi, "' is not numeric.")
    
    
    x[[i]] <- xi
  }
  
  
  # Convert from indices to values (now that ALL prefixes are removed).
  if (!is.null(within))  attr(x, 'within')  <- within  <- unique(x[within])
  if (!is.null(between)) attr(x, 'between') <- between <- unique(x[between])
  if (!is.null(revsort)) attr(x, 'revsort') <- revsort <- unique(x[revsort])
  
  
  assign(var, x, pos = env)
  
  return (invisible(NULL))
}




validate_bool <- function (var, env = parent.frame(), evar = var, max = 1, null_ok = FALSE, default = NULL) {
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  if (is.null(x)) x <- default
  
  if (!is_logical(x) || anyNA(x))
    stop ("`", evar, "` must be TRUE or FALSE.")
  
  validate_var_length(evar, x, max, null_ok)
  assign(var, x, pos = env)
  
  return (invisible(NULL))
}




validate_model <- function (var = "model", env = parent.frame()) {
  
  model <- get(var, pos = env, inherits = FALSE)
  
  
  #________________________________________________________
  # Predefined regression models
  #________________________________________________________
  if (is_scalar_character(model))
    model <- switch(
      EXPR = match.arg(model, c('lm', 'log', 'gam')),
      lm   = list("stats::lm", list(formula = y ~ x)),
      log  = list("stats::lm", list(formula = y ~ log(x))),
      gam  = list("mgcv::gam", list(formula = y ~ s(x, bs = "cs"), method = "REML" )) )
  
  
  
  #________________________________________________________
  # Sanity check model function.
  #________________________________________________________
  stopifnot(is_list(model))
  stopifnot(eq(length(model), 2L))
  fun  <- model[[1]]
  args <- model[[2]]
  
  stopifnot(is_list(args))
  stopifnot(is_formula(args[['formula']]))
  stopifnot(is_null(args[['data']]))
  
  if (is.character(fun) && !is_na(fun))
    fun <- local({
      fn  <- fun
      fun <- strsplit(fun, '::', fixed = TRUE)[[1]]
      fun <- if (length(fun) == 1) get(fun) else getFromNamespace(fun[[2]], fun[[1]])
      stopifnot(is_function(fun))
      stopifnot(all(c('formula', 'data') %in% formalArgs(fun)))
      attr(fun, 'fn') <- fn
      return (fun)
    })
  stopifnot(is.function(fun))
  
  fn <- attr(fun, 'fn', exact = TRUE)
  stopifnot(is_scalar_character(fn) && !is_na(fn))
  
  
  assign(var, list(fun = fun, args = args), pos = env)
  
  return (invisible(NULL))
}




validate_var_choices <- function (
    var, choices, env = parent.frame(), evar = var, max = 1, 
    null_ok = FALSE, na_ok = FALSE, all_option = NULL ) {
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  if (!missing(all_option) && eq(x, all_option)) {
    x <- choices
    
  } else {
    
    if (!is_character(x) && !is_na(x))
      stop(sprintf(
        fmt = "`%s` must be character, not %s",
        evar, 
        paste(collapse = ", ", class(x)) ))
    
    pm <- pmatch(tolower(x), tolower(choices))
    
    invalid_na <- which(is.na(pm))
    if (na_ok) invalid_na %<>% setdiff(which(is.na(x)))
    
    if (length(invalid_na) > 0)
      stop(sprintf(
        fmt = "Invalid `%s` option(s): %s\nChoices are: %s",
        evar, 
        paste(collapse = ", ", x[invalid_na]), 
        paste(collapse = ", ", choices) ))
    
    x <- choices[pm]
  }
  
  
  if (!na_ok && anyNA(x)) stop ("`", evar, "` cannot be NA.")
  if (length(x) > max)    stop ("`", evar, "` cannot be length ", length(x),".")
  
  assign(var, x, pos = env)
  return (invisible(NULL))
}


validate_var_range <- function (
    var, range, env = parent.frame(), evar = var, max = 1, 
    null_ok = FALSE, na_ok = FALSE ) {
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  if (!na_ok && anyNA(x))     stop ("`", evar, "` cannot be NA.")
  if (!null_ok && is.null(x)) stop ("`", evar, "` cannot be NULL.")
  if (!is.numeric(x))         stop ("`", evar, "` must be numeric.")
  if (length(x) > max)        stop ("`", evar, "` cannot be length ", length(x),".")
  
  if (any(x < min(range) | x > max(range)))
    stop ("`", evar, "` must be between ", min(range), " and ", max(range), ".")
  
  return (invisible(NULL))
}


validate_var_length <- function (evar, x, max, null_ok, na_ok = FALSE) {
  
  if (!na_ok && anyNA(x))     stop ("`", evar, "` cannot be NA.")
  if (!null_ok && is.null(x)) stop ("`", evar, "` cannot be NULL.")
  if (length(x) > max)        stop ("`", evar, "` cannot be length ", length(x),".")
  
  return (invisible(NULL))
}









#________________________________________________________
#' Convert all *.by params into long form.
#' 
#' @noRd
#' @keywords internal
#'
#' Examples values for color.by:
#' 
#' "Body Site"
#' c("Body Site" = "okabe")
#' list("Body Site")
#' list("Body Site" = "okabe")
#' list("Body Site" = list(values = "okabe", name = "Swab Location"))
#' list("Body Site" = list(values = c("red", "green", "blue"), name = "Swab Location"))
#' list("Body Site" = list(values = c(Saliva = "red", Stool = "green"), name = "Swab Location"))
#' 
#' 
#' c("Body Site", "Sex")
#' c("Body Site" = "okabe", "Sex")
#' c("Body Site" = "okabe", "Sex" = "muted")
#' list("Body Site", "Sex")
#' 
validate_meta_aes <- function (var, env = parent.frame(), null_ok = FALSE, ...) {
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  biom    <- get("biom", pos = env, inherits = FALSE)
  choices <- biom$fields
  is_num  <- sapply(USE.NAMES = TRUE, choices, function (i) is.numeric(biom$metadata[[i]]))
  
  
  #________________________________________________________
  # Allow ".all" option.
  #________________________________________________________
  if (eq(x, ".all")) x <- choices[-1]
  
  
  #________________________________________________________
  # Initial pass - convert `x` to a named list.
  #________________________________________________________
  results <- list()
  for (i in seq_along(x)) {
    
    key <- names(x)[[i]]
    val <- x[[i]]
    
    if (isTRUE(nzchar(key))) {
      results[[key]] <- val
    } else {
      stopifnot(is_scalar_character(val))
      results[[val]] <- list()
    }
  }
  
  
  #________________________________________________________
  # Remove prefixes and confirm metadata fields exist.
  #________________________________________________________
  col_names <- names(results)
  validate_meta('col_names', evar = var, null_ok = null_ok, ...)
  
  if (any(duplicated(names(results))))
    stop ("Duplicate metadata fields specified in `", var, "`.")
  
  names(results) <- col_names
  for (i in c("within", "between", "revsort"))
    attr(results, i) <- attr(col_names, i, exact = TRUE)
  
  remove("col_names", "i")
  
  
  
  #________________________________________________________
  # Check each metadata field's specs.
  #________________________________________________________
  for (col_name in names(results)) {
    
    col_spec <- results[[col_name]]
    col_type <- ifelse(is_num[[col_name]], "num", "cat")
    
    
    #________________________________________________________
    # Convert short-hand specification to long-form.
    #________________________________________________________
    result <- list()
    if (is_list(col_spec))           { result             <- col_spec
    } else if (is_palette(col_spec)) { result[['colors']] <- col_spec
    } else if (col_type == "num")    { 
      if (is.character(col_spec))    { result[['colors']] <- col_spec
      } else                         { result[['range']]  <- col_spec }
    } else                           { result[['values']] <- col_spec }
    
    
    #________________________________________________________
    # Only allow 'range' if numeric; 'values' if categorical.
    #________________________________________________________
    if (hasName(result, 'range')) {
      
      if (!eq(col_type, "num"))
        stop ("Can't apply `range` to non-numeric '", col_name, "' column.")
      
      if (!is.numeric(result[['range']]))
        stop ("Argument for '", col_name, "' `range` must be numeric.")
    }
    
    if (hasName(result, 'values')) {
      
      if (!eq(col_type, "cat"))
        stop ("Can't filter by `values` on non-categorical '", col_name, "' column.")
      
      if (!is.character(result[['values']]) && !eq(var, "shape.by"))
        stop ("Values for '", col_name, "' must be a character vector.")
    }
    
    
    attr(result, 'col_name') <- col_name
    attr(result, 'col_type') <- col_type
    
    results[[col_name]] <- result
  }
  
  
  if (length(results) == 0)
    results <- NULL
  
  
  assign(var, results, pos = env)
  return (invisible(NULL))
}



#________________________________________________________
#' Combines metadata variable's comparison attributes 
#' together with 'within' and 'between'.
#' 
#' Updates caller's `within` and `between` with attr(,'within') and 
#' attr(,'between') from variables named in vars. Ensures all within and 
#' between metadata column are categorical. Checks that within and between 
#' are not overlapping sets.
#' 
#' @noRd
#' @keywords internal
#' 

validate_meta_cmp <- function (vars, env = parent.frame()) {
  
  biom    <- get('biom',    pos = env, inherits = FALSE)
  within  <- get('within',  pos = env, inherits = FALSE)
  between <- get('between', pos = env, inherits = FALSE)
  
  for (i in mget(vars, env)) {
    within  %<>% c(attr(i, 'within',  exact = TRUE))
    between %<>% c(attr(i, 'between', exact = TRUE))
  }
  
  within  %<>% unique()
  between %<>% unique()
  
  validate_meta('within',  null_ok = TRUE, max = Inf, col_type = "cat")
  validate_meta('between', null_ok = TRUE, max = Inf, col_type = "cat")
  
  if (any(within %in% between))
    stop( "Metadata field name '", paste0(intersect(within, between), collapse = ", "), 
          "' cannot be set as both a within (==) and between (!=) grouping.")
  
  
  assign('within',  within,  pos = env)
  assign('between', between, pos = env)
  
  return (invisible(NULL))
}


