
cpu_count <- function () {
  
  #--------------------------------------------------------------
  # The number of cores we should use
  #--------------------------------------------------------------
  
  usercores <- try(as.numeric(getOption('rbiom.max.threads')), silent=TRUE)
  machcores <- try(parallel::detectCores(), silent=TRUE)
  consensus <- try(min(na.omit(c(usercores, machcores))), silent=TRUE)
  
  consensus <- if (!is(consensus, "try-error")) consensus else 1
  ncores    <- if (length(consensus) == 1)      consensus else 1
  ncores    <- if (is.numeric(ncores))          ncores    else 1
  ncores    <- if (ncores > 0 & ncores < Inf)   ncores    else 1
  ncores    <- as.integer(ncores)
  
  return (ncores)
}


configCluster <- function (nTasks=NA, pb=NULL) {
  
  ncores <- cpu_count()
  
  #--------------------------------------------------------------
  # Detect the currently available cluster
  #--------------------------------------------------------------
  
  cl <- list(
    'ready' = foreach::getDoParRegistered(),
    'name'  = foreach::getDoParName(),
    'cores' = foreach::getDoParWorkers()
  )
  
  
  #--------------------------------------------------------------
  # Start a new parallel backend
  #--------------------------------------------------------------
  
  if ( !identical(cl$ready, TRUE)    | 
       !identical(cl$name, "doSNOW") | 
       !identical(cl$cores, ncores) ) {
    
    if (identical(cl$ready, TRUE))
      cat(file=stderr(), "Replacing existing parallel backend.\n")
    
    msg <- sprintf("Setting up %i core cluster.\n", ncores)
    if (interactive()) cat(msg)
    try(pb$set(detail=msg), silent=TRUE)
    
    doSNOW::registerDoSNOW(parallel::makePSOCKcluster(ncores))
  }
  
  
  #--------------------------------------------------------------
  # Partition the tasks
  #--------------------------------------------------------------
  
  if (is.na(nTasks)) {
    
    res <- list(
      nTasks = ncores * 10, 
      ncores = ncores,
      sets   = seq_len(ncores * 10),
      nSets  = ncores * 10,
      opts   = NULL
    )
  } else {
  
    res <- list(
      nTasks = nTasks, 
      ncores = ncores,
      sets   = unlist(parallel::splitIndices(nTasks, min(nTasks, ncores * 10))),
      nSets  = min(nTasks, ncores * 10),
      opts   = NULL
    )
  }
  
  if (!is.null(pb)) {
    pb$set(0)
    res$opts <- list(progress = function (i) pb$set(i/res$nSets))
  }
  
  
  
  return (res)
    
}


progressBar <- function (progressbar=NULL, detail=NULL) {
  
  if (!is(progressbar, 'Progress')) {
    progressbar <- list(
      'set' = function (...) { invisible(NULL) },
      'inc' = function (...) { invisible(NULL) }
    )
  }
  
  if (!is.null(detail))
    try(progressbar$set(1, detail=detail), silent=TRUE)
  
  return (progressbar)
}




# withCluster <- function (expr, nTasks=NULL, env=NULL) {
#   
#   if (is.null(env)) env <- parent.frame()
#   
#   #--------------------------------------------------------------
#   # Autodetect and/or register a parallel backend
#   #--------------------------------------------------------------
#   
#   if (is.null(foreach:::.foreachGlobals[['data']])) {
#     
#     if (!is.numeric(ncores))
#       ncores <- min(getOption('rbiom.max.threads'), parallel::detectCores())
#       
#     if (!is.na(ncores) & ncores > 1) {
#       cl <- parallel::makeCluster(ncores)
#       doSNOW::registerDoSNOW(cl)
#     } else {
#       ncores <- 1
#       foreach::registerDoSEQ()
#     }
#     
#   } else {
#     ncores <- length(foreach:::.foreachGlobals[['data']])
#   }
#   
#   
#   
#   #--------------------------------------------------------------
#   # Partition the tasks
#   #--------------------------------------------------------------
#   if (!is.null(nTasks)) {
#     env$nTasks <- nTasks
#     env$sets   <- parallel::splitIndices(nTasks, min(nTasks, ncores * 20))
#     env$nSets  <- length(env$sets)
#     env$opts   <- list(progress = function (i) env$setProgress(value=i/env$nSets))
#   }
#   
#   
#   
#   #--------------------------------------------------------------
#   # Execute the code
#   #--------------------------------------------------------------
#   env$ncores <- ncores
#   expr       <- substitute(expr)
#   result     <- eval(expr, envir=env)
#   
#   
#   
#   #--------------------------------------------------------------
#   # Stop any cluster we've created, as is doesn't clean up 100%
#   #--------------------------------------------------------------
#   
#   # if (exists("cl", inherits=FALSE))
#   #   parallel::stopCluster(cl)
#   
#   
#   return (result)
#     
# }



# withProgress <- function (expr, progressbar=NULL, env=NULL) {
#   
#   if (is.null(env)) env <- parent.frame()
#   
#   expr <- substitute(expr)
#   
#   if (identical(progressbar, TRUE)) {
#     if (is(try(getDefaultReactiveDomain(), silent=TRUE), "ShinySession")) {
#       
#       env$setProgress <- shiny::setProgress
#       return (shiny::withProgress(expr, env=env))
#       
#     } else {
#       
#       env$setProgress <- function (value=NULL, message=NULL) {
#         cat(sprintf("%s\n", message))
#       }
#       
#       return (eval(expr, envir=env))
#     }
#   } else {
#     
#     env$setProgress <- function (value=NULL, message=NULL) {
#       invisible(NULL)
#     }
#     return (eval(expr, envir=env))
#   }
# }






