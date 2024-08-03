
stats_wilcox <- function (
    df, resp = attr(df, 'response'), stat.by = NULL, split.by = NULL, 
    level = 0.95, alt = "!=", mu = 0, p.adj = 'fdr' ) {
  
  
  regr <- NULL
  test <- "wilcox"
  fit  <- "lm"
  
  stats_validate()
  
  alternative <- switch(alt, '!=' = "two.sided", '>' = "greater", '<' = "less")
  
  # for CRAN check only
  estimate <- p.value <- conf.int <- statistic <- .stat.by <- NULL
  
  func     <- NULL
  pairwise <- NULL
  
  if (is.null(stat.by)) {
    
    pairwise <- FALSE # CRAN check requires consistent func signature
    func <- function (data, pair) {
      
      stats::wilcox.test(
          formula     = .resp ~ 1, 
          data        = data, 
          conf.int    = TRUE, 
          conf.level  = level, 
          exact       = FALSE, 
          alternative = alternative, 
          mu          = mu ) %>%
        with(tibble(
          row.names = NULL, 
          .mean     = estimate, 
          .h1       = factor(paste(alt, mu)), 
          .p.val    = p.value,
          .lower    = conf.int[[1]], 
          .upper    = conf.int[[2]], 
          .n        = nrow(data), 
          .stat     = statistic ))
      
    }
    
    code <- glue("
        stats <- wilcox.test(
            formula     = .resp ~ 1, 
            data        = data, 
            conf.int    = TRUE, 
            conf.level  = {level}, 
            exact       = FALSE, 
            alternative = '{alternative}', 
            mu          = {mu} ) %>%
          with(tibble(
            row.names = NULL, 
            .mean     = estimate, 
            .h1       = factor('{paste(alt, mu)}'), 
            .p.val    = p.value,
            .lower    = conf.int[[1]], 
            .upper    = conf.int[[2]], 
            .n        = nrow(data), 
            .stat     = statistic )) ")
    
    
  } else {
    
    pairwise <- TRUE # CRAN check requires consistent func signature
    func <- function (data, pair) {
        
        stats::wilcox.test(
            formula    = .resp ~ .stat.by, 
            data       = subset(data, .stat.by %in% pair), 
            conf.int   = TRUE, 
            conf.level = level, 
            exact      = FALSE ) %>%
          with(tibble(
            .mean.diff = estimate, 
            .h1        = factor('!= 0'), 
            .p.val     = p.value,
            .lower     = conf.int[[1]], 
            .upper     = conf.int[[2]], 
            .n         = sum(data$.stat.by %in% pair), 
            .stat      = statistic ))

    }
    code <- glue("
        stats <- wilcox.test(
            formula    = .resp ~ .stat.by, 
            data       = subset(data, .stat.by %in% pair), 
            conf.int   = TRUE, 
            conf.level = {level}, 
            exact      = FALSE ) %>%
          with(tibble(
            .mean.diff = estimate, 
            .h1        = factor('!= 0'), 
            .p.val     = p.value,
            .lower     = conf.int[[1]], 
            .upper     = conf.int[[2]], 
            .n         = sum(data$.stat.by %in% pair), 
            .stat      = statistic ))" )
    
  }
  
  
  stats <- stats_run(df, stat.by, split.by, func, pairwise) %>% aa(code = code)
  stats <- stats_finalize(stats, df, regr = NULL, resp, stat.by, split.by, fit = "lm", p.adj)
  
  attr(stats, 'tbl_sum')[['Model']] %<>% sub("^lm", "wilcox.test", .)
  
  return (stats)
}



stats_kruskal <- function (
    df, resp = attr(df, 'response'), stat.by, split.by = NULL, p.adj = 'fdr') {
  
  
  regr <- NULL
  test <- "kruskal"
  fit  <- "lm"
  
  
  stats_validate()
  
  if (is.null(stat.by))
    cli_abort("`stat.by` required when `test` = 'kruskal'.")
  
  
  pairwise <- FALSE
  func <- function (data, pair) {
    
    statistic <- p.value <- parameter <- NULL # for CRAN check only
    
    data %>% 
      stats::kruskal.test(.resp ~ .stat.by, .) %>%
      with(tibble(
        .stat  = statistic, 
        .h1    = factor('> 0'),
        .p.val = p.value, 
        .n     = nrow(data), 
        .df    = parameter ))
    
  }
  
  code <- glue("
    stats <- data %>% 
      stats::kruskal.test(.resp ~ .stat.by, .) %>%
      with(tibble(
        .stat  = statistic, 
        .h1    = factor('> 0'), 
        .p.val = p.value, 
        .n     = nrow(data), 
        .df    = parameter ))" )
  
  
  
  stats <- stats_run(df, stat.by, split.by, func, pairwise) %>% aa(code = code)
  stats <- stats_finalize(stats, df, regr = NULL, resp, stat.by, split.by, fit = "lm", p.adj, pairwise = FALSE)
  
  attr(stats, 'tbl_sum')[['Model']] %<>% sub("^lm", "kruskal.test", .)
  
  return (stats)
}
