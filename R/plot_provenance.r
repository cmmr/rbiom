

#________________________________________________________
# Log ggplot command history, complete with arguments.
# Called in .onLoad to assign function names in pkg env.
#________________________________________________________

ggwrap <- function (pkg, fn) {
  fun <- do.call(`::`, list(pkg, fn))
  
  assign(fn, pos = ENV, function (..., .indent = 0, .display = NULL) {
    res <- fun(...)
    
    if (is.null(.display))
      .display <- sprintf(
        fmt = "%s(%s)", 
        ifelse(pkg %in% c('ggplot2', 'grid'), fn, paste0(pkg, '::', fn)), 
        as.args(list(...), fun = fun, indent = .indent) )
    
    attr(res, 'display') <- .display
    
    return (res)
  })
}


#________________________________________________________
# Here, we'll call e.g. rep() using .rep() since we can't
# overwrite existing functions.
#________________________________________________________
# .rep <- function (..., .indent = 0) {
#   res  <- base::rep(...)
#   args <- as.args(list(...), fun = base::rep, indent = .indent)
#   attr(res, 'display') <- sprintf("rep(%s)", args)
#   return (res)
# }

basewrap <- function (pkg, fn) {
  assign(paste0(".", fn), pos = ENV, function (..., .indent = 0, .display = NULL) {
    fun  <- do.call(`::`, list(pkg, fn))
    res  <- fun(...)
    args <- as.args(list(...), fun = fun, indent = .indent)
    
    if (is.null(.display)) .display <- sprintf("%s(%s)", fn, args)
    attr(res, 'display') <- .display
    
    return (res)
  })
}

#________________________________________________________
# Add a layer to a list of layers.
#________________________________________________________

ggpush <- function (gglayers, gglayer) {
  
  if (is.null(attr(gglayer, 'display')))
    browser()
  
  gglayers[[length(gglayers) + 1]] <- gglayer
  return (gglayers)
}


#________________________________________________________
# Combine a list of logged commands into a plot.
#________________________________________________________

ggbuild <- function (gglayers) {
  
  p   <- NULL
  cmd <- NULL
  
  for (i in seq_along(gglayers)) {
    
    gglayer <- gglayers[[i]]
    
    # In case this layer was built by initLayer / setLayer
    if (!is.null(fun <- attr(gglayer, 'function', exact = TRUE)))
      gglayer <- do.call(fun, c(gglayer, '.indent' = 4))
    
    if (is.null(p)) {
      p   <- gglayer
      cmd <- attr(gglayer, 'display')
    } else {
      p   <- ggplot2::`%+%`(p, gglayer)
      cmd <- sprintf("%s +\n  %s", cmd, attr(gglayer, 'display'))
    }
  }
  
  attr(p, 'display') <- NULL
  attr(p, 'cmd')     <- cmd
  
  p$plot_env <- emptyenv()
  
  return (p)
}
