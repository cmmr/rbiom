

#________________________________________________________
# Log ggplot command history, complete with arguments.
# Called in .onLoad to assign function names in pkg env.
#________________________________________________________

cmd_wrap <- function (pkg, fn, env = NULL) {
  
  if (is_null(env)) env <- ENV
  
  if (nzchar(system.file(package = pkg))) {
  
    fun    <- getFromNamespace(x = fn, ns = pkg)
    pkgfn  <- ifelse(pkg %in% c('ggplot2', 'grid'), fn, paste0(pkg, '::', fn))
    
    newfun <- function (..., .indent = 0, .display = NULL, .lhs = NULL) {
      
      dots <- list(...)
      res  <- fun(...)
      
      
      if (is_null(.display)) {
        argstr   <- as.args(dots, fun = fun, indent = .indent)
        .display <- sprintf(fmt = "%s(%s)", pkgfn, argstr)
      }
      
      if (is_null(.lhs)) {
        attr(res, 'display') <- .display
      } else {
        attr(res, 'display') <- .lhs
        attr(res, 'code')    <- paste(collapse = "\n", c(
            if (length(dots) > 0) attr(dots[[1]], 'code', exact = TRUE) else NULL,
            sprintf("%s <- %s", .lhs, .display) )) %>% 
          add_class('rbiom_code')
      }
      
      
      return (res)
    }
    newfun %<>% aa(fn = fn, pkgfn = pkgfn, formalArgs = formalArgs(fun))
    assign(x = fn, value = newfun, pos = env)
    
    
  } else {
    
    # Case when an optional dependency is not installed.
    #________________________________________________________
    
    assign(fn, pos = env, function (...) {
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

basewrap <- function (pkg, fn, env = NULL) {
  
  if (is_null(env)) env <- ENV
  
  assign(paste0(".", fn), pos = env, function (..., .indent = 0, .display = NULL) {
    fun  <- do.call(`::`, list(pkg, fn))
    res  <- fun(...)
    args <- as.args(list(...), fun = fun, indent = .indent)
    
    if (is_null(.display)) .display <- sprintf("%s(%s)", fn, args)
    attr(res, 'display') <- .display
    
    return (res)
  })
}


