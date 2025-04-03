
FUNCS <- new.env()

P <- function (pkgfn) {
  
  if (env_has(FUNCS, pkgfn))
    return (FUNCS[[pkgfn]])
  
  
  parts <- strsplit(pkgfn, '::', fixed = TRUE)[[1]]
  pkg   <- parts[[1]]
  fn    <- parts[[2]]
  
  require_package(pkg, reason = paste0('to use `', fn, '()`'))
  
  fun <- getFromNamespace(x = fn, ns = pkg)
  fn  <- ifelse(pkg %in% c('base', 'ggplot2', 'grid'), fn, pkgfn)
  
  
  newfun <- function (..., .indent = 0) {
    
    dots <- list(...)
    res  <- fun(...)
    
    argstr  <- as.args(dots, fun = fun, indent = .indent)
    display <- sprintf(fmt = "%s(%s)", fn, argstr)
    
    attr(res, 'display') <- display
    
    return (res)
  }
  
  attr(newfun, 'display')    <- fn
  attr(newfun, 'formalArgs') <- formalArgs(fun)
  
  FUNCS[[pkgfn]] <- newfun
  
  return (newfun)
}


