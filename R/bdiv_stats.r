#' Test beta diversity vs categorical or numeric metadata.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_dist_test
#' @inherit documentation_default
#' @inherit documentation_stats_return return
#' 
#' @family beta_diversity
#' @family stats_tables
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     bdiv_stats(biom, stat.by = "Body Site", split.by = "Sex")
#'       
#'     bdiv_stats(biom, stat.by = "BMI", bdiv = c("bray", "unifrac"))
#'     
#'     # The R code used to compute the stats is in $code.
#'     tbl <- bdiv_stats(biom, stat.by = "Sex")
#'     tbl$code

bdiv_stats <- function (
    biom, stat.by, bdiv="Bray-Curtis", weighted = TRUE, tree = NULL, 
    test = "adonis2", within = NULL, between = NULL, 
    split.by = NULL, seed = 0, permutations = 999, p.adj = "fdr" ) {
  
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)
  
  params <- eval_envir(environment())
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity checks.
  #________________________________________________________
  with(params, {
    
    validate_bdiv(max = Inf)
    validate_bool('weighted', max = Inf)
    validate_meta('stat.by',  cmp = TRUE)
    validate_meta('split.by', cmp = TRUE, max = Inf, null_ok = TRUE)
    
    # Validates and appends to `within` and `between`.
    validate_meta_cmp(c('stat.by', 'split.by'))
    
    test  <-   match.arg(tolower(test), c("adonis2", "mrpp"))
    p.adj %<>% match.arg(p.adjust.methods)
    
    stopifnot(is_scalar_integerish(seed)         && !is_na(seed))
    stopifnot(is_scalar_integerish(permutations) && !is_na(permutations))
    
  })
  
  
  #________________________________________________________
  # Iteratively compute dm and permutational stats.
  #________________________________________________________
  stats <- with(params, bdply(
    biom   = biom, 
    vars   = split.by, 
    iters  = list(weighted = weighted, bdiv = bdiv),
    prefix = TRUE,
    FUN    = function (b, weighted, bdiv) {
      groups <- pull(b, stat.by)
      dm     <- bdiv_distmat(biom = b, bdiv = bdiv, weighted = weighted, tree = tree)
      distmat_stats(dm = dm, groups = groups, test = test, seed = seed, permutations = permutations)
    })) %>%
    transform(.adj.p = stats::p.adjust(.p.val, params$p.adj)) %>%
    plyr::arrange(.p.val) %>%
    as_rbiom_tbl()
  
  
  #________________________________________________________
  # Table header.
  #________________________________________________________
  attr(stats, 'tbl_sum') <- with(params, c(
    'Test' = glue::glue(switch(
      EXPR = test,
      adonis2 = "vegan::adonis2(~ {coan(stat.by)}, permutations = {permutations})",
      mrpp    = "vegan::mrpp(grouping = {coan(stat.by)}, permutations = {permutations})" ))))
  
  
  
  #________________________________________________________
  # Attach human-readable commands of varying complexity.
  #________________________________________________________
  attr(stats, 'code') <- with(params, {
    
    if (is.null(split.by) && length(weighted) == 1 && length(bdiv) == 1) {
      
      #________________________________________________________
      # Simple version without any iterations.
      #________________________________________________________
      glue::glue(
        .sep = "\n",
        "dm <- bdiv_distmat(biom = biom, bdiv = '{bdiv}', weighted = {weighted})",
        "grouping <- pull(biom, {as.args(list(stat.by))})[attr(dm, 'Labels')]",
        "",
        "set.seed({seed})",
        "",
        switch(
          EXPR = test,
          adonis2 = "vegan::adonis2(formula = dm ~ grouping, permutations = {permutations}) %>%",
          mrpp    = "vegan::mrpp(dat = dm, grouping = grouping, permutations = {permutations}) %>%" ),
        "  vegan::permustats() %>%", 
        "  summary() %>%", 
        "  with(data.frame(.stat = statistic, .z = z, .p.val = p)) %>%", 
        "  tryCatch(", 
        "    error   = function (e) data.frame(.stat=NA, .z=NA, .p.val=NA), ", 
        "    warning = function (w) data.frame(.stat=NA, .z=NA, .p.val=NA) ) %>%", 
        "  data.frame(row.names = NULL, .n = attr(dm, 'Size'), .) %>%",
        "  data.frame({as.args(list(.weighted=weighted, .bdiv=bdiv))}, .) %>%",
        "  transform(.adj.p = stats::p.adjust(.p.val, '{p.adj}')) %>%",
        "  plyr::arrange(.p.val)" )
      
      
    } else {
      
      #________________________________________________________
      # With iterations.
      #________________________________________________________
      glue::glue(
        .sep = "\n",
        "bdply(",
        "  biom   = biom, ",
        "  vars   = {as.args(list(split.by))}, ",
        "  iters  = list({as.args(list(weighted = weighted, bdiv = bdiv))}), ",
        "  prefix = TRUE, ",
        "  FUN    = function (b, weighted, bdiv) {{",
        "",
        "    dm <- bdiv_distmat(biom = b, bdiv = bdiv, weighted = weighted)",
        "    grouping <- pull(b, stat.by)[attr(dm, 'Labels')]",
        "",
        "    set.seed({seed})",
        "",
        switch(
          EXPR = test,
          adonis2 = "    vegan::adonis2(formula = dm ~ grouping, permutations = {permutations}) %>%",
          mrpp    = "    vegan::mrpp(dat = dm, grouping = grouping, permutations = {permutations}) %>%" ),
        "      vegan::permustats() %>%", 
        "      summary() %>%", 
        "      with(data.frame(.stat = statistic, .z = z, .p.val = p)) %>%", 
        "      tryCatch(", 
        "        error   = function (e) data.frame(.stat=NA, .z=NA, .p.val=NA), ", 
        "        warning = function (w) data.frame(.stat=NA, .z=NA, .p.val=NA) ) %>%", 
        "      data.frame(row.names = NULL, .n = attr(dm, 'Size'), .)",
        "",
        "}}) %>%",
        "  transform(.adj.p = stats::p.adjust(.p.val, '{p.adj}')) %>%",
        "  plyr::arrange(.p.val)" )
      
    }
    
  }) %>% add_class('rbiom_code')
  
  
  
  set_cache_value(cache_file, stats)
  
  return (stats)
}
