

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
  if (inherits(x, 'phylo'))  return (invisible(NULL))
  
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
    n <- length(choices)
    x <- ifelse(
      test = x >= 0, 
      yes  = choices[pmin(x + 1,     n)], 
      no   = choices[pmax(x + 1 + n, 1)] )
  }
  
  validate_var_choices('x', choices, null_ok = null_ok, evar = var, ...)
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



validate_layers <- function (var = "layers", choices, env = parent.frame()) {
  
  x <- get(var, pos = env, inherits = FALSE)
  
  if (anyNA(x))         cli_abort("`{var}` cannot have any NA values.")
  if (!is_character(x)) cli_abort("`{var}` must be character, not {.type {x}}.")
  
  if (is.null(names(choices))) names(choices) <- substr(choices, 1, 1)
  
  layer_names <- NULL
  for (i in tolower(x[nchar(x) > 0]))
    layer_names <- c(layer_names, local({
      
      if (i %in% choices)        return (i)
      if (i %in% names(choices)) return (choices[[i]])
      
      pm <- pmatch(i, choices)
      ii <- strsplit(i, '')[[1]]
      
      if (length(x) > 1 && !is.na(pm)) return (choices[[pm]])
      if (all(ii %in% names(choices))) return (choices[ii])
      if (!is.na(pm))                  return (choices[[pm]])
      return (choices[intersect(ii, names(choices))])
    }))
  
  layer_names <- layer_names[!is.na(layer_names)]
  
  
  if (length(layer_names) == 0)
    cli_abort("Invalid `{var}` argument: {.val {x}}.")
  
  assign(var, layer_names, pos = env)
}



validate_string <- function (
    var, env = parent.frame(), evar = var, 
    null_ok = FALSE, na_ok = FALSE, 
    if_true = TRUE, if_false = FALSE, if_na = NA, if_null = NULL ) {
  
  x <- get(var, pos = env, inherits = FALSE)
  
  if (is_empty(x)) x <- NULL
  if (is_true(x))  x <- if_true
  if (is_false(x)) x <- if_false
  if (is_null(x))  x <- if_null
  if (is_na(x))    x <- if_na
  
  if (length(x) > 1)          stop ("`", evar, "` can't be longer than 1 value.")
  if (is_null(x) && !null_ok) stop ("`", evar, "` cannot be NULL.")
  if (anyNA(x)   && !na_ok)   stop ("`", evar, "` cannot be NA.")
  
  assign(var, x, pos = env)
  return (invisible(NULL))
}


validate_var_choices <- function (
    var, choices, env = parent.frame(), evar = var, max = 1, 
    null_ok = FALSE, na_ok = FALSE, dne_ok = FALSE, all_option = NULL ) {
  
  if (!hasName(env, var)) {
    if (dne_ok) return (invisible(NULL))
    cli_abort("{.var {evar}} does not exist.")
  }
  
  x <- get(var, pos = env, inherits = FALSE)
  if (length(x) == 0) x <- NULL # Convert character(0) to NULL
  
  if (is.null(x) && null_ok) {
    assign(var, x, pos = env)
    return (invisible(NULL))
  }
  
  if (is.data.frame(choices)) choices <- colnames(choices)
  
  if (!missing(all_option) && eq(x, all_option)) {
    x <- choices
    
  } else {
    
    if (!is_character(x) && !is_na(x))
      cli_abort("{.var {evar}} must be character, not {.type {x}}.")
    
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
    var, range = c(-Inf, Inf), env = parent.frame(), evar = var, 
    n = NA, int = FALSE, null_ok = FALSE, na_ok = FALSE, dne_ok = FALSE ) {
  
  if (!hasName(env, var)) {
    if (dne_ok) return (invisible(NULL))
    cli_abort("{.var {evar}} does not exist.")
  }
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  if (!na_ok && anyNA(x))     stop ("`", evar, "` cannot be NA.")
  if (!null_ok && is.null(x)) stop ("`", evar, "` cannot be NULL.")
  
  if (!is.numeric(x) && !all(is.na(x)))
    cli_abort("{.var {evar}} must be numeric, not {.type {x}}.")
  
  if (!is.na(n) && length(x) != n)
    stop ("`", evar, "` must be length ", n, ", not ", length(x), ".")
    
  if (isTRUE(int) && any(!is.na(x) & x %% 1 != 0))
    stop ("`", evar, "` must be whole numbers.")
  
  if (any(!is.na(x) & (x < min(range) | x > max(range))))
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
#' Combines metadata variable's comparison attributes 
#' together with `within` and `between`.
#' 
#' * Updates `within` and `between` in caller's environment. 
#' * Variables named in `vars` are stripped of '!=' and '==' prefixes.
#' 
#' @noRd
#' @keywords internal
#' 

validate_var_cmp <- function (vars, env = parent.frame()) {
  
  within  <- get('within',  pos = env, inherits = FALSE)
  between <- get('between', pos = env, inherits = FALSE)
  
  
  for (var in vars) {
    
    x <- get(var, pos = env, inherits = FALSE)
    if (!is.character(x)) next
    
    for (i in seq_along(x)) {
      xi <- x[[i]]
      if (is.na(xi)) next
      if (startsWith(xi, '==')) within  %<>% c(substr(xi, 3, nchar(xi)))
      if (startsWith(xi, '!=')) between %<>% c(substr(xi, 3, nchar(xi)))
      if (startsWith(xi, '==') || startsWith(xi, '!=')) {
        x[[i]] <- substr(xi, 3, nchar(xi))
      }
    }
    
    assign(var, x, pos = env)
  }
  
  
  assign('within',  within,  pos = env)
  assign('between', between, pos = env)
}




#________________________________________________________
#' Ensures that `var` corresponds to a metadata field.
#' 
#' @noRd
#' @keywords internal
#' 

validate_biom_field <- function (
    var, env = parent.frame(), evar = var, 
    null_ok = FALSE, max = 1, col_type = NULL ) {
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  biom <- get('biom', pos = env, inherits = FALSE)
  df   <- env$biom$metadata
  
  validate_df_field(
    var      = 'x', 
    evar     = var, 
    null_ok  = null_ok, 
    max      = max, 
    col_type = col_type, 
    force    = FALSE, 
    drop_na  = FALSE )
  
  assign(var, x, pos = env)
  
  return (invisible(NULL))
}


validate_df_field <- function (
    var, env = parent.frame(), evar = var, 
    null_ok = FALSE, max = 1, col_type = NULL, 
    force = TRUE, drop_na = TRUE ) {
  
  x <- get(var, pos = env, inherits = FALSE)
  if (is.null(x) && null_ok) return (invisible(NULL))
  
  df      <- get('df', pos = env, inherits = FALSE)
  choices <- colnames(df)
  is_num  <- if (!is.null(col_type)) lapply(df, is.numeric)
  
  
  if (eq(x, '.all')) {
    x <- setdiff(choices, '.sample')
    
    # Exit early if '.all' resolved to NULL
    if (length(x) == 0 && null_ok) {
      assign(var, NULL, pos = env)
      return (invisible(NULL))
    }
  }
  
  
  if (is_null(x))       cli_abort("`{evar}` cannot be NULL.")
  if (!is_character(x)) cli_abort("`{evar}` should be character, not {.type {x}}.")
  if (length(x) > max)  cli_abort("`{evar}` cannot be length ", length(x), ".")
  if (anyNA(x))         cli_abort("`{evar}` cannot be NA.")
  x <- trimws(x)
  if (!all(nzchar(x)))  cli_abort("`{evar}` cannot be ''.")
  
  
  for (i in seq_along(x)) {
    
    xi <- x[[i]]
    
    
    # Check for valid name and (optionally) type.
    if (!xi %in% choices)
      validate_var_choices('xi', choices, evar = evar)
    
    if (eq(col_type, "cat") && isTRUE(force)) {
      df[[xi]] %<>% refactor()
      is_num[[xi]] <- FALSE
    }
    
    if (eq(col_type, "cat") && is_num[[xi]])
      cli_abort("`{evar}` field '{xi}' is not categorical.")
    
    if (eq(col_type, "num") && !is_num[[xi]])
      cli_abort("`{evar}` field '{xi}' is not numeric.")
    
    if (isTRUE(drop_na))
      df <- df[!is.na(df[[xi]]),,drop=FALSE]
    
    x[[i]] <- xi
  }
  
  
  if (isTRUE(force) || isTRUE(drop_na))
    assign('df', df, pos = env)
  
  assign(var, unique(x), pos = env)
  
  return (invisible(NULL))
}


is_nz_string <- function (str) {
  isTRUE(is_scalar_character(str) && !is_na(str) && nzchar(trimws(str)))
}

is_url <- function (str) {
  if (!is_nz_string(str)) return (FALSE)
  isTRUE(grepl('^(http|ftp)s{0,1}://', tolower(str)))
}
