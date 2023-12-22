
#' Run non-parametric statistics on a data.frame.
#' 
#' @inherit documentation_test.ifelse
#' @inherit documentation_model.lm
#' @inherit documentation_default
#' @inherit documentation_stats_return return
#' 
#' @family stats_tables
#' 
#' @param df   A data.frame with columns named by \code{stat.by}, \code{resp}, 
#'        \code{regr}, and \code{split.by}. Required.
#' 
#' @param resp   The response (independent) numeric variable, such as taxa 
#'        abundance or alpha diversity.  Default: \code{attr(df, 'response')}
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     df <- taxa_table(biom, rank = "Family")
#'     stats_table(df, stat.by = "Body Site")
#'     
#'     df <- adiv_table(biom)
#'     stats_table(df, stat.by = "Sex", test = "pw_means")
#'     
#'     stats_table(df, stat.by = "Sex", regr = "BMI")
#' 

stats_table <- function (
    df, stat.by, resp = attr(df, 'response'), regr = NULL, 
    test = ifelse(is.null(regr), "means", "trends"), model = "lm", 
    level = 0.95, split.by = NULL, p.adj = 'fdr' ) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment())
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Sanity checks.
  #________________________________________________________
  validate_var_range('level', c(0,1))
  validate_var_choices('p.adj', choices = p.adjust.methods)
  
  stopifnot(is.data.frame(df) && nrow(df) > 0)
  choices <- colnames(df)
  validate_var_choices('stat.by',  choices, null_ok = TRUE)
  validate_var_choices('resp',     choices)
  validate_var_choices('regr',     choices, null_ok = TRUE)
  validate_var_choices('split.by', choices, null_ok = TRUE, max = Inf)
  remove("choices")
  
  test <- match.arg(
    arg     = trimws(tolower(test)),
    choices = if (is.null(regr)) {
        c("means", "pw_means")
      } else {
        c("predict", "terms", "fit", "means", "trends", 
          "es_means", "es_trends", "pw_means", "pw_trends" ) })
  
  if (is.null(stat.by) && !test %in% c('predict', 'fit', 'terms', 'means', 'trends'))
    stop ("Test '", test, "' requires a `stat.by` group.")
  
  
  
  #________________________________________________________
  # Response (`resp`) variable must be numeric.
  #________________________________________________________
  if (!is.numeric(df[[resp]]))
    stop("Column '", resp, "' must be numeric, not ", head(class(df[[resp]])), ".")
  
  df <- df[is.finite(df[[resp]]),,drop=FALSE] # Drop NA, Inf, NaN, etc
  
  
  
  #________________________________________________________
  # Coerce `stat.by` and `split.by` variables to factors.
  #________________________________________________________
  for (i in c(stat.by, split.by))
    if (!is.factor(df[[i]])) {
      if (!is.character(df[[i]]))
        cli_warn("Numeric column '{i}' is being used for categorical grouping.")
      df[[i]] %<>% as.factor()
    }
  remove(list = intersect("i", ls()))
  
  
  
  #________________________________________________________
  # Calculate boxplot-like or corrplot-like stats.
  #________________________________________________________
  if (is.null(regr)) {
    args  <- fun_params(stats_table_cat, as.list(environment()))
    stats <- do.call(stats_table_cat, args)
    
  } else {
    args  <- fun_params(stats_table_num, as.list(environment()))
    stats <- do.call(stats_table_num, args)
  }
  
  if (hasName(stats, '.p.val'))
    stats[['.adj.p']] <- p.adjust(stats[['.p.val']], p.adj)
  
  
  
  #________________________________________________________
  # Add in "data <- ..."
  #________________________________________________________
  if (!is.null(attr(df, 'cmd', exact = TRUE)))
    attr(stats, 'code') %<>% sprintf("data  <- %s\n%s", attr(df, 'cmd'), .)
  
  
  attr(stats, 'stats_table_args') <- args[names(args) != 'df']
  attr(stats, 'code') %<>% add_class("rbiom_code")
  
  
  #________________________________________________________
  # Enable accessing attributes with `$`.
  #________________________________________________________
  stats %<>% as_rbiom_tbl()
  
  
  set_cache_value(cache_file, stats)
  
  return (stats)
}
