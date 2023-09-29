

#________________________________________________________
# Log ggplot command history, complete with arguments.
# Called in .onLoad to assign function names in pkg env.
#________________________________________________________

cmd_wrap <- function (pkg, fn) {
  
  if (nzchar(system.file(package = pkg))) {
  
    fun    <- getFromNamespace(x = fn, ns = pkg)
    pkgfn  <- ifelse(pkg %in% c('ggplot2', 'grid'), fn, paste0(pkg, '::', fn))
    
    newfun <- function (..., .indent = 0, .display = NULL, .lhs = NULL) {
      
      res <- fun(...)
      
      if (is_null(.display)) {
        argstr   <- as.args(list(...), fun = fun, indent = .indent)
        .display <- sprintf(fmt = "%s(%s)", pkgfn, argstr)
      }
      
      if (is_null(.lhs)) {
        attr(res, 'display') <- .display
      } else {
        attr(res, 'cmd')     <- sprintf("%s <- %s", .lhs, .display)
        attr(res, 'display') <- .lhs
      }
      
      
      return (res)
    }
    newfun %<>% aa(fn = fn, pkgfn = pkgfn)
    assign(x = fn, value = newfun, pos = ENV)
    
    
  } else {
    
    # Case when an optional dependency is not installed.
    #________________________________________________________
    
    assign(fn, pos = ENV, function (...) {
      stop("R package '", pkg, "' must be installed to use this feature.")
    })
  }
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
    
    if (is_null(.display)) .display <- sprintf("%s(%s)", fn, args)
    attr(res, 'display') <- .display
    
    return (res)
  })
}


#________________________________________________________
# Add a layer to a list of layers.
#________________________________________________________

ggpush <- function (gglayers, gglayer) {
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
    if (!is_null(fun <- attr(gglayer, 'function', exact = TRUE)))
      gglayer <- do.call(fun, c(gglayer, '.indent' = 4))
    
    if (is_null(p)) {
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
