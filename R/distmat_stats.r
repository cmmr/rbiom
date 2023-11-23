
#' Run statistics on a distance matrix vs a categorical or numeric variable.
#' 
#' @inherit documentation_dist_test
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family stats_tables
#'        
#' @return A data.frame with summary statistics from [vegan::permustats()]. 
#'         The columns are:
#'        \itemize{
#'          \item{\emph{.n} - }{ The size of the distance matrix. }
#'          \item{\emph{.stat} - }{ 
#'                The observed statistic. For mrpp, this is the overall 
#'                weighted mean of group mean distances. }
#'          \item{\emph{.z} - }{ 
#'                The difference of observed statistic and mean of permutations 
#'                divided by the standard deviation of permutations (also known 
#'                as z-values). Evaluated from permuted values without observed 
#'                statistic. }
#'          \item{\emph{.p.val} - }{ Probability calculated by \code{test}. }
#'        }\cr
#'         R commands for reproducing the results are in \code{$code} and 
#'         object history in \code{$history}.
#'         
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- sample_select(hmp50, 1:10)
#'     dm   <- bdiv_distmat(biom, 'unifrac')
#'     
#'     distmat_stats(dm, groups = sample_metadata(biom, 'Body Site'))
#'     
#'     distmat_stats(dm, groups = sample_metadata(biom, 'Age'))
#'     
#'     # See the R code used to calculate these statistics:
#'     stats <- distmat_stats(dm, groups = sample_metadata(biom, 'Age'))
#'     stats$code
#'

distmat_stats <- function (dm, groups, test = "adonis2", seed = 0, permutations = 999) {
  
  
  params  <- eval_envir(environment())
  history <- append_history('stats ', params)
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  #________________________________________________________
  # Sanity checking
  #________________________________________________________
  with(params, {
    test <- match.arg(tolower(test), c("adonis2", "mrpp"))
    stopifnot(is_scalar_integerish(seed) && !is.na(seed))
    stopifnot(is_scalar_integerish(permutations) && !is.na(permutations))
    stopifnot(!is.null(names(groups)))
    stopifnot(is(dm, 'dist'))
  })
  
  
  #________________________________________________________
  # Calculate distance matrix statistics.
  #________________________________________________________
  stats <- with(params, {
    
    grouping <- groups[attr(dm, 'Labels')]
    set.seed(seed)
    
    test %>%
      switch(
        adonis2 = vegan::adonis2(formula = dm ~ grouping, permutations = permutations),
        mrpp    = vegan::mrpp(dat = dm, grouping = grouping, permutations = permutations) ) %>%
      vegan::permustats() %>%
      summary() %>%
      with(data.frame(.stat = statistic, .z = z, .p.val = p)) %>%
      tryCatch(
        error   = function (e) data.frame(.stat=NA, .z=NA, .p.val=NA), 
        warning = function (w) data.frame(.stat=NA, .z=NA, .p.val=NA) ) %>%
      data.frame(row.names = NULL, .n = attr(dm, 'Size'), .) %>%
      as_rbiom_tbl()
  }) 
  
  
  #________________________________________________________
  # Present the code used to the user.
  #________________________________________________________
  attr(stats, 'code') <- with(params, glue::glue(
      .sep = "\n",
      "grouping <- groups[attr(dm, 'Labels')]",
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
      "  data.frame(row.names = NULL, .n = attr(dm, 'Size'), .)"
      
    )) %>% add_class('rbiom_code')
  
  
  #________________________________________________________
  # Return the statistics table.
  #________________________________________________________
  set_cache_value(cache_file, stats)
  return (stats)
  
}

