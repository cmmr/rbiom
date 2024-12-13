

###  Stats Helpers  -----------------
#____________________________________
#____________________________________


stats_validate <- function () {
  
  with(parent.frame(), {
    
    df <- df[,c(regr, resp, stat.by, split.by),drop = FALSE]
    df <- df[stats::complete.cases(df),,drop = FALSE]
    if (nrow(df) == 0) cli_abort("No rows after removing NA values.")
    
    
    #________________________________________________________
    # Rename columns to '.regr', '.resp', '.stat.by'
    #________________________________________________________
    
    if (!is.null(regr))
      names(df)[which(names(df) == regr)] <- ".regr"
    
    names(df)[which(names(df) == resp)] <- ".resp"
    
    if (!is.null(stat.by))
      names(df)[which(names(df) == stat.by)] <- ".stat.by"
    
  })
}




stats_run <- function (df, stat.by, split.by, func, pairwise) {
  
  
  if (isFALSE(pairwise)) {
    
    plyr::ddply(df, ply_cols(split.by), function (data) {
      
      tryCatch(
        error = function (e) data.frame()[1,], 
        expr  = suppressWarnings({
          
          func(data, NULL)
          
        }))
      
    }) %>%
      as_tibble()
    
    
  } else {
    
    plyr::ddply(df, ply_cols(split.by), function (data) {
      
      tryCatch(
        error = function (e) data.frame()[1,], 
        expr  = suppressWarnings({
          
          pairs <- combn(intersect(levels(data$.stat.by), data$.stat.by), m = 2)
          colnames(pairs) <- plyr::aaply(pairs, 2L, paste, collapse = ' - ')
          
          plyr::adply(pairs, 2L, .id = stat.by, function (pair) {
            
            tryCatch(
              error = function (e) data.frame()[1,], 
              expr  = suppressWarnings({
                
                func(data, pair)
                
              }))
            
          })
        }))
      
    }) %>%
      as_tibble()
    
  }
  
  
}




stats_finalize <- function (stats, df, regr, resp, stat.by, split.by, fit, p.adj, pairwise = !is.null(stat.by)) {
  
  code <- attr(stats, 'code', exact = TRUE)
  
  #________________________________________________________
  # Pairwise STAT.BY wrapper.
  #________________________________________________________
  if (!is.null(stat.by) && isTRUE(pairwise)) {
    code <- sub('stats <- ', '', code, fixed = TRUE)
    code <- paste0(
      "pairs <- combn(intersect(levels(data$.stat.by), data$.stat.by), m = 2)\n",
      "colnames(pairs) <- plyr::aaply(pairs, 2L, paste, collapse = ' - ')\n\n",
      "stats <- plyr::adply(pairs, 2L, .id = ", single_quote(stat.by), ", function (pair) {\n",
      "  tryCatch(error = function (e) data.frame()[1,], suppressWarnings({\n\n",
      "    ", gsub("\n", "\n    ", code), "\n\n",
      "  }))\n",
      "})", 
      ifelse(is.null(split.by), " %>% \n  as_tibble()", "")
    )
  }
  
  
  #________________________________________________________
  # SPLIT.BY wrapper.
  #________________________________________________________
  if (!is.null(split.by)) {
    code <- sub('stats <- ', '', code, fixed = TRUE)
    code <- paste0(
      "stats <- plyr::ddply(data, .(", paste(collapse = ', ', coan(split.by)), "), function (data) {\n",
      "  tryCatch(error = function (e) data.frame()[1,], suppressWarnings({\n\n",
      "    ", gsub("\n", "\n    ", code), "\n\n",
      "  }))\n",
      "}) %>% \n",
      "  tibble::as_tibble()"
    )
  }
  
  
  #________________________________________________________
  # Drop auto-generated column when split.by is NULL.
  #________________________________________________________
  if (hasName(stats, '.id'))
    stats <- dplyr::select(stats, -any_of('.id'))
  
  
  #________________________________________________________
  # The `regr` field should be first.
  #________________________________________________________
  if (!is.null(regr) && hasName(stats, regr) && names(stats)[[1]] != regr) {
    stats <- dplyr::relocate(stats, any_of(regr), .before = 1)
    code %<>% paste0(" %>% \n", glue("  dplyr::relocate({single_quote(regr)}, .before = 1)"))
  }
  
  
  #________________________________________________________
  # Correct for multiple comparisons and arrange by p.val.
  #________________________________________________________
  if (hasName(stats, '.p.val')) {
    
    .p.val <- NULL # for CRAN check only
    
    stats <- stats %>% 
      dplyr::mutate(.adj.p = p.adjust(.p.val, p.adj), .after = .p.val) %>% 
      dplyr::arrange(.p.val)
    
    code %<>% paste0(
      " %>% \n",
      "  dplyr::mutate(.adj.p = p.adjust(.p.val, '", p.adj, "'), .after = .p.val) %>% \n", 
      "  dplyr::arrange(.p.val)" )
  }
  
  
  
  #________________________________________________________
  # Prepend the renaming of data fields.
  #________________________________________________________
  code <- paste0(
    "data %<>% dplyr::rename(", 
    as.args(indent = 2, as.list(c(
      .regr = regr,
      .resp = resp,
      .stat.by = stat.by
    ))), 
    ")\n\n", code )
  
  
  #________________________________________________________
  # Prepend code for creating the original data.
  #________________________________________________________
  if (!is.null(attr(df, 'code', exact = TRUE))) {
    code <- paste0(attr(df, 'code', exact = TRUE), "\n\n", code)
    
  } else if (!is.null(attr(df, 'cmd', exact = TRUE))) {
    code <- paste0("data <- ", attr(df, 'cmd', exact = TRUE), "\n\n", code)
  }
  
  
  
  #________________________________________________________
  # Enable special printing.
  #________________________________________________________
  stats <- as_rbiom_tbl(stats)
  attr(stats, 'tbl_sum') <- c('Model' = stats_model_str(fit, stat.by, regr, resp))
  attr(stats, 'code')    <- code
  
  
  return (stats)
}




stats_formula <- function (fit, stat.by = FALSE, regr = TRUE, resp = ".resp", k = NULL) {
  
  stopifnot(fit %in% c('lm', 'log', 'gam'))
  
  resp <- coan(resp)
  
  if (isTRUE(stat.by))   stat.by <- ".stat.by"
  if (isFALSE(stat.by))  stat.by <- NULL
  if (!is.null(stat.by)) stat.by <- coan(stat.by)
  
  if (isTRUE(regr))   regr <- ".regr"
  if (isFALSE(regr))  regr <- NULL
  if (!is.null(regr)) regr <- coan(regr)
  
  k <- ifelse(is.null(k), '', paste(", k =", k))
  
  
  frm <- if (is.null(regr) && is.null(stat.by)) {
    
    "{resp} ~ 1"
    
  } else if (is.null(regr)) {
    
    "{resp} ~ {stat.by}"
  
  } else if (is.null(stat.by)) {
    
    switch(
      EXPR = fit,
      'lm'  = "{resp} ~ {regr}", 
      'log' = "{resp} ~ log({regr})", 
      'gam' = "{resp} ~ s({regr}, bs='cs'{k})" )
    
  } else {
    
    switch(
      EXPR = fit,
      'lm'  = "{resp} ~ {regr} * {stat.by}", 
      'log' = "{resp} ~ log({regr}) * {stat.by}", 
      'gam' = "{resp} ~ s({regr}, by={stat.by}, bs='cs'{k}) + {stat.by}" )
  }
  
  as.formula(glue(frm), env = baseenv())
}



stats_fit_model <- function (df, fit, stat.by = FALSE, regr = ".regr", resp = ".resp") {
  
  fun  <- switch(fit, lm = stats::lm, log = stats::lm, gam = mgcv::gam )
  k    <- floor(max(1, min(10, nrow(df) / 5)))
  frm  <- stats_formula(fit, stat.by, regr, resp, k = k)
  args <- if (fit == "gam") list(method = "REML") else list()
  
  do.call(fun, c(list(frm, data = df), args))
  
}



stats_model_str <- function (fit, stat.by = FALSE, regr = ".regr", resp = ".resp") {
  
  stopifnot(fit %in% c('lm', 'log', 'gam'))
  
  frm  <- stats_formula(fit, stat.by, regr, resp)
  args <- if (eq(fit, "gam")) list(method = "REML") else list()
  if (eq(fit, "log")) fit <- "lm"
  
  glue("{fit}({as.args(c(list(frm), args))})")
}


stats_model_code <- function (var, fit, stat.by) {
  
  if (fit == "gam") {
    
    frm <- stats_formula(fit, stat.by, k = 'gam_k')
    
    glue("
    
      gam_k <- floor(max(1, min(10, nrow({var}) / 5)))
      gam_f <- {format(frm)}
      model <- mgcv::gam(formula = gam_f, data = {var}, method = 'REML')")
    
  } else {
    
    frm <- stats_formula(fit, stat.by)
    glue("model <- stats::lm({format(frm)}, {var})")
  }
  
  
}


