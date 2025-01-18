
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
#'        \describe{
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
#'        }
#'         R commands for reproducing the results are in \code{$code}.
#'         
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     hmp10        <- hmp50$clone()
#'     hmp10$counts <- hmp10$counts[,1:10]
#'     
#'     dm <- bdiv_distmat(hmp10, 'unifrac')
#'     
#'     distmat_stats(dm, groups = pull(hmp10, 'Body Site'))
#'     
#'     distmat_stats(dm, groups = pull(hmp10, 'Age'))
#'     
#'     # See the R code used to calculate these statistics:
#'     stats <- distmat_stats(dm, groups = pull(hmp10, 'Age'))
#'     stats$code
#'

distmat_stats <- function (dm, groups, test = "adonis2", seed = 0, permutations = 999) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env()
  cache_file <- get_cache_file('distmat_stats', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity checking
  #________________________________________________________
  params <- list2env(params)
  with(params, {
    test <- match.arg(tolower(test), c("adonis2", "mrpp", "none"))
    validate_seed()
    stopifnot(is_scalar_integerish(permutations) && !is.na(permutations))
    stopifnot(!is.null(names(groups)))
    stopifnot(inherits(dm, 'dist'))
  })
  if (params$test == "none") return (NULL)
  
  
  #________________________________________________________
  # Calculate distance matrix statistics.
  #________________________________________________________
  stats <- with(params, {
    
    grouping <- groups[attr(dm, 'Labels')]
    
    
    # Preserve current .Random.seed
    oldseed <- if (exists(".Random.seed")) .Random.seed else NULL
    set.seed(seed)
    
    result <- test %>%
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
    
    if (!is.null(oldseed)) .Random.seed <- oldseed
    
    return (result)
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
  attr(stats, 'cmd') <- current_cmd('distmat_stats')
  set_cache_value(cache_file, stats)
  return (stats)
  
}

