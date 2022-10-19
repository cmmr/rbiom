

#________________________________________________________
# Log ggplot command history, complete with arguments.
# Called in .onLoad to assign function names in pkg env.
#________________________________________________________

ggwrap <- function (pkg, fn) {
  assign(fn, pos = ENV, function (..., .indent = 0, .display = NULL) {
    fun  <- do.call(`::`, list(pkg, fn))
    res  <- fun(...)
    args <- as.args(list(...), fun = fun, indent = .indent)
    if (!pkg %in% c('ggplot2', 'grid')) fn <- paste0(pkg, '::', fn)
    
    if (is.null(.display)) .display <- sprintf("%s(%s)", fn, args)
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
  
  res <- NULL
  cmd <- NULL
  
  for (i in seq_along(gglayers)) {
    
    gglayer <- gglayers[[i]]
    
    if (is.null(res)) {
      res <- gglayer
      cmd <- attr(gglayer, 'display')
    } else {
      res <- ggplot2::`%+%`(res, gglayer)
      cmd <- sprintf("%s +\n  %s", cmd, attr(gglayer, 'display'))
    }
  }
  attr(res, 'cmd') <- cmd
  return (res)
}
