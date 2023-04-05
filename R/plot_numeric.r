
# See also ggpmisc package
# https://docs.r4photobiology.info/ggpmisc/articles/model-based-annotations.html


plot_numeric <- function (biom, x, y, layers = "pr", mode, params) {
  
  #-----------------------------------------------
  # Convert user's `layers` spec to layer names.
  #-----------------------------------------------
  layer_names <- local({
    layerlist <- c( 'p' = "point", 'r' = "regression")
    
    layer_match(layers, choices = layerlist, default = "pr") %>%
      c('ggplot', ., 'xaxis', 'yaxis', 'labs', 'theme', 'theme_bw')
  })
  
  
  # #-----------------------------------------------
  # # Default arguments for the layers
  # #-----------------------------------------------
  # layers[['labs']][['y']] <- as.vector(y)
  # 
  # 
  # 
  # 
  # #-----------------------------------------------
  # # The x, y, etc to plot
  # #-----------------------------------------------
  # for (i in unique(c(color.by, shape.by, facet.by)))
  #   biom[['metadata']][[i]] %<>% as.factor()
  # 
  # md     <- unique(c(x, color.by, shape.by, facet.by))
  # ggdata <- distill(biom = biom, metric = y, md = md, safe = TRUE)
  # names(ggdata) <- sub(attr(ggdata, 'response', exact = TRUE), ".y", names(ggdata))
  # 
  # 
  # #-----------------------------------------------
  # # All the layers will share a single data frame
  # #-----------------------------------------------
  # ggdata[['.src']] <- "points"
  # 
  # 
  # #-----------------------------------------------
  # # Trend Line and Statistics
  # #-----------------------------------------------
  # if ("regression" %in% names(layers)) {
  #   
  #   
  #   #-----------------------------------------------
  #   # Parse the formula from 'model'
  #   #-----------------------------------------------
  #   tryCatch({
  #     formula <- if (!is(model, 'formula')) as.formula(model) else model
  #     stopifnot(is(formula, 'formula'))
  #   }, error = function (e) stop("Invalid value for `formula`. ", e) )
  #   
  #   
  #   #-----------------------------------------------
  #   # Does 'regr' looks like a regression function?
  #   #-----------------------------------------------
  #   tryCatch({
  #     stopifnot(is(regr, 'function'))
  #     if (!all(c('formula', 'data') %in% formalArgs(regr)))
  #       stop("function does not have parameters 'formula' and 'data'.")
  #   }, error = function (e) stop("Invalid value for `regr`. ", e) )
  #   
  #   
  #   
  #   #-----------------------------------------------
  #   # Parameters for stat_smooth()
  #   #-----------------------------------------------
  #   layers[['regression']][['formula']] <- formula
  #   layers[['regression']][['method']]  <- regr
  #   
  #   if (is_null(ci)) {
  #     layers[['regression']][['se']] <- FALSE
  #     
  #   } else if (is.numeric(ci)) {
  #     if (ci < 75 || ci > 100)
  #       stop("`ci` must be between 75 and 100. Value given = ", ci)
  #     layers[['regression']][['se']]    <- TRUE
  #     layers[['regression']][['level']] <-  ci / 100
  #     
  #   } else {
  #     stop("`ci` must be either NULL or a number between 75 and 100. Value given = ", ci)
  #   }
  #   
  #   
  #   #-----------------------------------------------
  #   # Run Statistics
  #   #-----------------------------------------------
  #   
  #   # change `y ~ x` into `Shannon ~ BMI`
  #   formula <- do.call(substitute, list(expr = model, env = list(x=as.name(x), y=as.name(y))))
  #   
  #   ply_by <- unique(intersect(names(ggdata), c(color.by, shape.by, facet.by)))
  #   ply_df <- ggdata[,intersect(names(ggdata), c(".value", all.vars(formula), ply_by)),drop=F]
  #   names(ply_df) <- sub(".value", y, names(ply_df))
  #   
  #   stats_df <- plyr::ddply(ply_df, backtick(ply_by), function(x) {
  #     tryCatch(
  #       error   = function (e) return (NULL),
  #       warning = function (w) return (NULL),
  #       {
  #         regr(formula, data = x) %>% broom::glance() %>% signif(3)
  #       })
  #   })
  #   if ('p.value' %in% names(stats_df))
  #     stats_df[['adj.p.value']] <- p.adjust(stats_df[['p.value']], method = p.adj)
  #   
  #   
  #   #-----------------------------------------------
  #   # Add Annotations
  #   #-----------------------------------------------
  #   anno <- layer_match(anno, c(names(stats_df), "method"), "mpr")
  #   
  #   if ('method' %in% anno)
  #     layers[['labs']][['caption']] <- paste0("Regression formula: ", fun_toString(regr), "(", format(formula), ")" )
  #   
  #   
  #   anno <- setdiff(anno, "method")
  #   if (isTRUE(length(anno) > 0 && nrow(stats_df) > 0)) {
  #     stats_df[['.label']] <- apply(stats_df[,anno,drop=F], 1L, function (x) paste(names(x), sep=" = ", x, collapse="; "))
  #     
  #     if (is_null(facet.by)) {
  #       layers[['labs']][['subtitle']] <- stats_df[['.label']]
  #       
  #     } else {
  #       
  #       #-----------------------------------------------
  #       # Add stats to overall data frame
  #       #-----------------------------------------------
  #       layers[['stats_text']] <- list(
  #         data    = ~ subset(.x, .src == "stats"),
  #         mapping = list(label = ".label"),
  #         x       = -Inf, 
  #         y       =  Inf,
  #         hjust   = -0.04
  #       )
  #       
  #       
  #       #-----------------------------------------------
  #       # df will be appended to ggdata
  #       #-----------------------------------------------
  #       df <- stats_df[,c(ply_by, '.label')]
  #       df[['.src']] <- "stats"
  #       
  #       
  #       #-----------------------------------------------
  #       # Prevent labels from colliding with each other
  #       #-----------------------------------------------
  #       if (!is_null(color.by) && !identical(color.by, facet.by)) {
  #         ypos <- as.numeric(df[[color.by]])
  #         df[['.vjust']] <- 1.5 * ypos
  #         layers[['stats_text']][['mapping']] %<>% c(list(vjust = ".vjust"))
  #         layers[['yaxis']][['expand']] <- structure(
  #           expansion(mult = c(0, 0.08) * max(ypos)), 
  #           'display' = sprintf("expansion(mult = c(0, %.2f))", 0.08 * max(ypos)) )
  #         
  #       } else {
  #         layers[['stats_text']][['vjust']] <- 1.5
  #         layers[['yaxis']][['expand']] <- structure(
  #           expansion(mult = c(0, 0.08)), 
  #           'display' = sprintf("expansion(mult = c(0, 0.08))") )
  #       }
  #       
  #       ggdata %<>% append_df(df)
  #     }
  #     
  #     stats_df <- stats_df[,setdiff(names(stats_df), '.label'),drop=F]
  #   }
  # }
  # 
  # 
  # #--------------------------------------------------------------
  # # Set up the global aesthetics
  # #--------------------------------------------------------------
  # aes_args <- list(x = as.name(x), y = ".value")
  # 
  # 
  # 
  # #--------------------------------------------------------------
  # # One-off tweaks
  # #--------------------------------------------------------------
  # layers[['labs']][['y']] %<>% sub(
  #   fixed       = TRUE,
  #   pattern     = "OTUs Diversity", 
  #   replacement = "Observed OTUs" )
  
  
  
  p <- layers %>%
    plot_numeric_facets() %>%
    plot_numeric_stats() %>%
    plot_build()
  
  return (p)
}
