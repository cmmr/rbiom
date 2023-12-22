
#' Run Kruskal-Wallis or Mann-Whitney statistics.
#' 
#' All arguments should be pre-validated by [stats_table()].
#' 
#' @noRd
#' @keywords internal
#' 

stats_table_cat <- function (df, stat.by, resp, test, level, split.by) {
  
  stopifnot(is_string(test, c("means", "pw_means")))
  
  
  #________________________________________________________
  # Construct the formula for the stats method.
  #________________________________________________________
  frm <- as.formula(paste(backtick(resp), "~", backtick(stat.by)))
  
  
  #________________________________________________________
  # Pairwise stats: Wilcox Rank Sum, aka Mann-Whitney.
  #________________________________________________________
  if (test == "pw_means") {
    
    pairs <- combn(levels(df[[stat.by]]), m = 2, simplify = FALSE)
    names(pairs) <- sapply(pairs, paste, collapse = " - ")
    
    result <- plyr::ddply(df, ply_cols(split.by), function (d) {
      plyr::ldply(pairs, .id = stat.by, function (pair) {
        
        pair <- d[d[[stat.by]] %in% pair,,drop=FALSE]
        
        wilcox.test(
            formula    = frm, 
            data       = pair, 
            conf.int   = TRUE, 
            conf.level = level ) %>%
          with(data.frame(
            row.names = NULL,
            .n        = nrow(pair),
            .stat     = statistic, 
            .estimate = estimate, 
            .lower    = conf.int[[1]], 
            .upper    = conf.int[[2]], 
            .p.val    = p.value )) %>%
          signif(digits = 3) %>%
          suppressWarnings() %>%
          tryCatch(error = function (e) data.frame()[1,])
        
      })
    }) %>% as_tibble()
    
    attr(result, 'pairs') <- pairs
    
    
    if (is.null(split.by) || length(split.by) == 0) {
      
      attr(result, 'code') <- glue::glue("
        pairs <- combn(levels(data${coan(stat.by)}), m = 2, simplify = FALSE)
        names(pairs) <- sapply(pairs, paste, collapse = ' - ')
        
        stats <- plyr::ldply(pairs, .id = {single_quote(stat.by)}, function (pair) {{
          
          pair <- data[data${coan(stat.by)} %in% pair,,drop=FALSE]
          
          wilcox.test({format(frm)}, data = pair, conf.int = TRUE, conf.level = {level}) %>%
            with(data.frame(
              row.names = NULL,
              .n        = nrow(pair),
              .stat     = statistic, 
              .estimate = estimate, 
              .lower    = conf.int[[1]], 
              .upper    = conf.int[[2]], 
              .p.val    = p.value )) %>%
            suppressWarnings() %>%
            tryCatch(error = function (e) data.frame()[1,])
        }})"
      )
      
    } else {
      
      attr(result, 'code') <- glue::glue("
        pairs <- combn(levels(data${coan(stat.by)}), m = 2, simplify = FALSE)
        names(pairs) <- sapply(pairs, paste, collapse = ' - ')
        
        stats <- plyr::ddply(data, single_quote(split.by), function (d) {{
          plyr::ldply(pairs, .id = {single_quote(stat.by)}, function (pair) {{
            
            pair <- d[d${coan(stat.by)} %in% pair,,drop=FALSE]
            
            wilcox.test({format(frm)}, data = pair, conf.int = TRUE, conf.level = {level}) %>%
              with(data.frame(
                row.names = NULL,
                .n        = nrow(pair),
                .stat     = statistic, 
                .estimate = estimate, 
                .lower    = conf.int[[1]], 
                .upper    = conf.int[[2]], 
                .p.val    = p.value )) %>%
              suppressWarnings() %>%
              tryCatch(error = function (e) data.frame()[1,])
          }})
        }})"
      )
      
    }
  }
  
  
  
  #________________________________________________________
  # Groupwise stats: Kruskal-Wallis
  #________________________________________________________
  if (test == "means") {
    
    result <- plyr::ddply(df, ply_cols(split.by), function (d) {
      
      stats::kruskal.test(formula = frm, data = d) %>%
        with(data.frame(
          row.names = NULL, 
          .n     = nrow(d), 
          .stat  = statistic, 
          .df    = parameter, 
          .p.val = p.value )) %>%
        tryCatch(error = function (e) data.frame()[1,])
      
    }) %>% as_tibble()
    
    
    if (is.null(split.by) || length(split.by) == 0) {
      
      attr(result, 'code') <- glue::glue("
        stats <- stats::kruskal.test({format(frm)}, data = data) %>%
          with(data.frame(
           .n     = nrow(data), 
           .stat  = statistic, 
           .df    = parameter, 
           .p.val = p.value ))"
      )
      
    } else {
      
      attr(result, 'code') <- glue::glue("
        stats <- plyr::ddply(data, {as.args(list(split.by))}, function (d) {{
          
          stats::kruskal.test({format(frm)}, data = d) %>%
              
            with(data.frame(
              row.names = NULL, 
              .n     = nrow(d), 
              .stat  = statistic, 
              .df    = parameter, 
              .p.val = p.value )) %>%
              
            tryCatch(error = function (e) data.frame()[1,])
          
        }})"
      )
      
    }
    
  }
  
  
  # When split.by is NULL, '.id' = NA appears.
  result <- result[,colnames(result) != ".id"]
  
  
  class(result) %<>% c('rbiom_tbl', .)
  attr(result, 'tbl_sum') <- c('Test' = switch(
    EXPR = test, 
    means    = paste0("kruskal.test(", format(frm), ")."), 
    pw_means = paste0("wilcox.test(",  format(frm), ").") ))
  
  return (result)
}
