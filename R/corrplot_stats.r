

#________________________________________________________
# Computes p-values for categorical differences.
#________________________________________________________
corrplot_stats <- function (layers) {
  
  params <- attr(layers, 'params', exact = TRUE)
  xcol   <- attr(layers, 'xcol',   exact = TRUE)
  ycol   <- attr(layers, 'ycol',   exact = TRUE)
  ggdata <- attr(layers, 'data',   exact = TRUE)
  
  color.by  <- names(params[['color.by']])
  facet.by  <- params[['facet.by']]
  ci        <- params[['ci']]
  model     <- params[['model']]
  p.adj     <- params[['p.adj']]
  level     <- ifelse(is.numeric(ci), ci / 100, 0.95)
  
  color.by   <- unique(setdiff(color.by, facet.by))
  facet.by   <- unique(facet.by)
  emm_specs  <- if (!is.null(color.by)) color.by else facet.by
  emm_by     <- if (!is.null(color.by)) facet.by else NULL
  predictors <- c(emm_specs, emm_by)
  
  stopifnot(is_string(p.adj, stats::p.adjust.methods))
  
  
  #________________________________________________________
  # If anything goes wrong, skip stats.
  #________________________________________________________
  stats <- tryCatch(
    error = function (e) { warning(e); return (NULL) }, 
    expr  = local({
      
      if (is_null(model))                       return (NULL)
      if (identical(params[['stats']], "none")) return (NULL)
      
      
      #________________________________________________________
      # Convert from geom_smooth formula to manual stats.
      #________________________________________________________
      model_function <- model[[1]]
      model_args     <- model[[2]]
      
      model_args[['formula']] <- local({
        replacements <- list(x = as.symbol(xcol), y = as.symbol(ycol))
        f <- eval(do.call(substitute, list(model_args[['formula']], replacements)))
        if (is_null(color.by)) return (update(f, ~ . + 0))
        update(f, sprintf("~ . * %s + 0", backtick(color.by)))
      })
      
      err <- function (...) data.frame()[1,]
      
      
      
      results <- plyr::dlply(ggdata, ply_cols(facet.by), function (df) {
        
        
        #________________________________________________________
        # Can't run stats on just a single factor level.
        #________________________________________________________
        if (!is_null(color.by) && length(unique(df[[color.by]])) < 2) {
          
          x <- data.frame()[1,]
          if (length(facet.by) > 0)
            for (i in facet.by) x[[i]] <- df[1,i]
            
          return (list(
            'df' = df, 'fit' = x, 'terms' = x, 
            'emmeans'  = x, 'emm_pairs' = x, 'emm_eff_size' = x, 
            'emtrends' = x, 'emt_pairs' = x, 'emt_eff_size' = x ))
        }
        
        
        res <- list()
        
        #________________________________________________________
        # Create a model just for this subset (facet) of ggdata.
        #________________________________________________________
        model_args[['data']] <- df
        m <- do.call(model_function, model_args, envir = baseenv())
        
        
        
        #________________________________________________________
        # Stats that come into play for color.by
        #________________________________________________________
        if (length(color.by) > 0) {
          
          
          #________________________________________________________
          # Look at estimated means.
          #________________________________________________________
          emm <- try(silent = TRUE, emmeans::emmeans(object = m, specs = c(xcol, color.by), level = level, infer = TRUE))
          res[['emmeans']]      <- tryCatch(summary(object = emm, adjust = 'none'),                                            error = err, warning = err)
          res[['emm_pairs']]    <- tryCatch(pairs(x = emm, adjust = 'none', simple = color.by),                                error = err, warning = err)
          res[['emm_eff_size']] <- tryCatch(eff_size(object = emm, sigma = sigma(m), edf = df.residual(m), simple = color.by), error = err, warning = err)
          
          
          #________________________________________________________
          # Look at linear trends (slopes).
          #________________________________________________________
          emt <- try(silent = TRUE, emmeans::emtrends(object = m, specs = color.by, var = xcol, level = level, infer = TRUE))
          res[['emtrends']]     <- tryCatch(summary(object = emt, adjust = 'none'),                         error = err, warning = err)
          res[['emt_pairs']]    <- tryCatch(pairs(x = emt, adjust = 'none'),                                error = err, warning = err)
          res[['emt_eff_size']] <- tryCatch(eff_size(object = emt, sigma = sigma(m), edf = df.residual(m)), error = err, warning = err)
          
          for (i in head(which(endsWith(names(res[['emtrends']]), ".trend")), 1))
            names(res[['emtrends']])[[i]] <- "slope"
          
        }
        
        
        
        #________________________________________________________
        # Basic stats for all models.
        #________________________________________________________
        res[['df' ]]   <- broom::augment(m, data = df)
        res[['terms']] <- broom::tidy(x = m)
        res[['fit'  ]] <- broom::glance(x = m)
        
        
        #________________________________________________________
        # Prepend facet columns to the data frames; drop xcol.
        #________________________________________________________
        res <- lapply(res, function (x) {
          x <- as.data.frame(x)
          for (i in facet.by) x[[i]] <- df[1,i]
          x <- x[,setdiff(unique(c(facet.by, colnames(x))), xcol),drop=FALSE]
          return (x)
        })
        
        
        return (res)
      })
      
      
      #________________________________________________________
      # Combine all per-facet tables into one.
      #________________________________________________________
      stats <- list()
      for (i in unique(as.vector(sapply(results, names)))) {
        stats[[i]] <- plyr::rbind.fill(lapply(results, `[[`, i))
        
        if (hasName(stats[[i]], "p.value"))
          stats[[i]][['adj.p']] <- p.adjust(stats[[i]][['p.value']], p.adj)
        
        for (j in seq_len(ncol(stats[[i]])))
          if (is.numeric(stats[[i]][[j]]))
            stats[[i]][[j]] <- signif(stats[[i]][[j]], 3)
      }
      
      
      
      #________________________________________________________
      # Attach stats R commands as attr(,'cmd') attributes.
      #________________________________________________________
      
      stats <- local({
        
        model_fn      <- attr(model_function, 'fn', exact = TRUE)
        model_arg_str <- as.args(model_args)
        facet_str     <- as.args(list(facet.by))
        color.by_str  <- glue::double_quote(color.by)
        xcol_str      <- glue::double_quote(xcol)
        
        emm_template <- ifelse(
          test = is.null(facet.by), 
          yes  = paste0(
            sprintf("model <- %s(%s, data = data)\n", model_fn, model_arg_str),
            sprintf("emm   <- emmeans(object = model, specs = c(%s, %s), level = %s, infer = TRUE)\n", xcol_str, color.by_str, level),
            sprintf("stats <- {cmd}") ),
          no   = paste0(
            sprintf("stats <- plyr::ddply(data, %s, function (df) {\n", facet_str),
            sprintf("  model <- %s(%s, data = df)\n", model_fn, model_arg_str),
            sprintf("  emm   <- emmeans(object = model, specs = c(%s, %s), level = %s, infer = TRUE)\n", xcol_str, color.by_str, level),
            sprintf("  {cmd}\n}") ))
        
        if (hasName(stats, 'emmeans'))
          attr(stats[['emmeans']],  'cmd') <- emm_template %>%
            sub("{cmd}", "summary(object = emm, adjust = 'none')", ., fixed = TRUE) %>%
            paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
        
        if (hasName(stats, 'emm_pairs'))
          attr(stats[['emm_pairs']],  'cmd') <- emm_template %>%
            sub("{cmd}", sprintf("pairs(x = emm, simple = %s, adjust = 'none')", color.by_str), ., fixed = TRUE) %>%
            paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
        
        if (hasName(stats, 'emm_eff_size'))
          attr(stats[['emm_eff_size']], 'cmd') <- emm_template %>%
            sub("{cmd}", "eff_size(object = emm, sigma = sigma(model), edf = df.residual(model))", ., fixed = TRUE)
        
        
        
        emt_template <- ifelse(
          test = is.null(facet.by), 
          yes  = paste0(
            sprintf("model <- %s(%s, data = data)\n", model_fn, model_arg_str),
            sprintf("emt   <- emtrends(object = model, specs = %s, var = %s, level = %s, infer = TRUE)\n", color.by_str, xcol_str, level),
            sprintf("stats <- {cmd}") ),
          no   = paste0(
            sprintf("stats <- plyr::ddply(data, %s, function (df) {\n", facet_str),
            sprintf("  model <- %s(%s, data = df)\n", model_fn, model_arg_str),
            sprintf("  emt   <- emtrends(object = model, specs = %s, var = %s, level = %s, infer = TRUE)\n", color.by_str, xcol_str, level),
            sprintf("  {cmd}\n}") ))
        
        if (hasName(stats, 'emtrends'))
          attr(stats[['emtrends']],  'cmd') <- emt_template %>%
            sub("{cmd}", "summary(object = emt, adjust = 'none')", ., fixed = TRUE) %>%
            paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
        
        if (hasName(stats, 'emt_pairs'))
          attr(stats[['emt_pairs']],  'cmd') <- emt_template %>%
            sub("{cmd}", "pairs(x = emt, adjust = 'none')", ., fixed = TRUE) %>%
            paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
        
        if (hasName(stats, 'emt_eff_size'))
          attr(stats[['emt_eff_size']], 'cmd') <- emt_template %>%
            sub("{cmd}", "eff_size(object = emt, sigma = sigma(model), edf = df.residual(model))", ., fixed = TRUE)
        
        
        broom_template <- ifelse(
          test = is.null(facet.by), 
          yes  = paste0(
            sprintf("model <- %s(%s, data = data)\n", model_fn, model_arg_str),
            sprintf("stats <- {cmd}") ),
          no   = paste0(
            sprintf("stats <- plyr::ddply(data, %s, function (df) {\n", facet_str),
            sprintf("  model <- %s(%s, data = df)\n", model_fn, model_arg_str),
            sprintf("  {cmd}\n}") ))
        
        if (hasName(stats, 'df'))
          attr(stats[['df']],  'cmd') <- broom_template %>%
            sub("{cmd}", "broom::augment(x = model, data = {data})", ., fixed = TRUE) %>%
            sub("{data}", ifelse(is.null(facet.by), "data", "df"), ., fixed = TRUE)
        
        if (hasName(stats, 'terms'))
          attr(stats[['terms']],  'cmd') <- broom_template %>%
            sub("{cmd}", "broom::tidy(x = model)", ., fixed = TRUE) %>%
            paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
        
        if (hasName(stats, 'fit'))
          attr(stats[['fit']], 'cmd') <- broom_template %>%
            sub("{cmd}", "broom::glance(x = model)", ., fixed = TRUE) %>%
            paste0("\nstats[['adj.p']] <- p.adjust(stats[['p.value']], '", p.adj, "')")
        
        
        return (stats)
      })
      
      
      #________________________________________________________
      # Noteworthy arguments used here. Only `unbox`-ables.
      #________________________________________________________
      attr(stats, 'args') <- list(
        xcol     = xcol, 
        ycol     = ycol, 
        color.by = color.by )
      
      return (stats)
    }))
  
  
  #________________________________________________________
  # Use p.top to retain only most significant taxa.
  #________________________________________________________
  if (isTRUE(is.finite(params[['p.top']])))
    if (!is.null(stats[[params[['stats']]]][['p.value']])) {
      
      taxa <- if (params[['p.top']] >= 1) {
        rank(plyr::daply(stats[[params[['stats']]]], '.taxa', function (s) {
          min(c(Inf, s[['p.value']]), na.rm = TRUE) }), ties.method = "first")
        
      } else {
        plyr::daply(stats[[params[['stats']]]], '.taxa', function (s) {
          min(c(Inf, s[['adj.p']]), na.rm = TRUE) })
      }
      
      taxa   <- unique(names(taxa[taxa <= params[['p.top']]]))
      ggdata <- ggdata[ggdata[['.taxa']] %in% taxa,,drop=FALSE]
      ggdata[['.taxa']] %<>% factor()
      
      for (i in names(stats)) {
        df <- stats[[i]]
        if (!hasName(df, ".taxa")) next
        df <- df[df[['.taxa']] %in% taxa,,drop=FALSE]
        df[['.taxa']] %<>% factor()
        stats[[i]] <- df
      }
      
      for (i in names(attributes(ggdata))) {
        df <- attr(ggdata, i, exact = TRUE)
        if (!is.data.frame(df))    next
        if (!hasName(df, ".taxa")) next
        df <- df[df[['.taxa']] %in% taxa,,drop=FALSE]
        df[['.taxa']] %<>% factor()
        attr(ggdata, i) <- df
      }
      
      remove("taxa")
    }
  
  
  
  # if (!is.null(stats[['fit']])) {
  #   
  #   stats_text <- sprintf(
  #     fmt = "*p* = %s; *R<sup>2</sup>* = %s; *F* = %s",
  #     format(stats[['fit']][['p.value']],   digits=3), 
  #     format(stats[['fit']][['r.squared']], digits=3), 
  #     format(stats[['fit']][['statistic']], digits=3) )
  #   
  #   
  # 
  #   #________________________________________________________
  #   # Add caption describing the model/formula.
  #   #________________________________________________________
  #   model_cmd <- local({
  #     
  #     fun  <- model[[1]]
  #     args <- model[[2]]
  #     
  #     fm <- capture.output(args[['formula']])[[1]]
  #     for (i in predictors)
  #       fm %<>% paste(sep = " + ", capture.output(as.symbol(i)))
  #     args[['formula']] <- structure(fm, display = fm)
  #     
  #     str <- sprintf("%s(%s)", attr(fun, "fn", exact = TRUE), as.args(args, fun = fun))
  #     
  #     # Ensure that nothing in the formula is interpreted as markdown syntax.
  #     entities <- c(
  #       setNames(paste0("&#", 33:42,    ";"), strsplit("!\"#$%&'()*", "")[[1]]),
  #       setNames(paste0("&#", c(60,62), ";"), strsplit("<>", "")[[1]]),
  #       setNames(paste0("&#", 91:96,    ";"), strsplit("[\\]^_`", "")[[1]]),
  #       setNames(paste0("&#", 123:126,  ";"), strsplit("{|}~", "")[[1]]) )
  #     for (i in seq_along(entities))
  #       str <- gsub(str, pattern = names(entities)[[i]], replacement = entities[[i]], fixed = TRUE)
  #     
  #     return (str)
  #   })
  #   
  #   methods_text <- ifelse(
  #     test = isFALSE(ci),
  #     yes  = sprintf("Curve fitted using %s", model_cmd),
  #     no   = sprintf("Curve and %g%% CI fitted using %s", ci, model_cmd) )
  #   
  #   subtitle <- sprintf("%s<br><span style='font-size:9pt'>%s</span>", stats_text[[1]], methods_text)
  #   setLayer("labs",  subtitle = subtitle)
  #   setLayer("theme", plot.subtitle = element_markdown(size = 11, lineheight = 1.2))
  #   
  # }
  
  
  params[['stats']] <- NULL
  params[['model']] <- NULL
  params[['p.adj']] <- NULL
  params[['p.top']] <- NULL
  
  attr(layers, 'data')   <- ggdata
  attr(layers, 'stats')  <- stats
  attr(layers, 'params') <- params
  
  
  return (layers)
}
