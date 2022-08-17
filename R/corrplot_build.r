

corrplot_build <- function (layers) {
  
  params <- attr(layers, 'params', exact = TRUE)
  
  
  #________________________________________________________
  # aes() mapping for color/fill
  #________________________________________________________
  if (!is.null(params[['color.by']])) {
    if (hasName(layers, 'point'))
      setLayer("point", "mapping|color" = params[['color.by']])
    
    if (hasName(layers, 'smooth')) {
      setLayer("smooth", "mapping|color" = params[['color.by']])
      if (!isFALSE(params[['ci']]))
        setLayer("smooth", "mapping|fill" = params[['color.by']])
    }
  }
  
  
  #________________________________________________________
  # Fade out the point when a curve is fitted.
  #________________________________________________________
  if (hasLayer("smooth") && hasLayer("point"))
    setLayer("point", size = 0.2, alpha = 0.5)
  
  
  #________________________________________________________
  # Define the curve and confidence interval.
  #________________________________________________________
  if (hasLayer("smooth")) {
    
    setLayer("smooth", params[c('method', 'formula')])
    params[grep("\\.(method|formula)$", names(params))] <- NULL
    
    ci <- params[['ci']]
    if (is.logical(ci)) setLayer("smooth", se    = ci)
    if (is.numeric(ci)) setLayer("smooth", level = ci / 100)
    params[grep("(^|\\.)(se|level)$", names(params))] <- NULL
  }
  
  
  #________________________________________________________
  # Add caption describing the model/formula.
  #________________________________________________________
  if (hasLayer("smooth")) {
    
    setLayer("labs", caption = local({
      
      caption <- ifelse(
        test = isFALSE(params[['ci']]), 
        yes  = "Curve fitted using", 
        no   = sprintf("Curve and %g%% CI fitted using", params[['ci']]) )
      
      model <- as.args(list(params[['method']]))
      model <- gsub('(\\w+\\:\\:|")', '', model) # remove quotes, pkg::
      
      if (identical(model, 'loess')) {
        caption %<>% paste("local polynomial regression")
        
      } else {
        if (identical(model, 'lm'))  model <- "linear model"
        if (identical(model, 'glm')) model <- "generalized linear model"
        if (identical(model, 'nlm')) model <- "non-linear model"
        if (identical(model, 'gam')) model <- "generalized additive model"
        
        caption %<>% paste(as.args(list(params[['formula']])), model)
      }
      
      return (paste0(caption, "."))
    }))
    
    setLayer("theme", plot.caption = element_text(face = "italic"))
  }
  
  
  
  
  attr(layers, 'params') <- params
  
  
  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  p <- layers %>% 
    plot_facets() %>% 
    plot_build()
  
  
  return (p)
}
