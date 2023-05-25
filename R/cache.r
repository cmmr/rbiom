
#' Speeds up repetitive computations by storing results of functions calls. 
#' 
#' @name with_cache
#' 
#' @param expr  An expression to cache the output of.
#' 
#' @return The result of \code{expr}.
#' @examples
#'     library(rbiom)
#'     
#'     system.time(x <- adiv_table(hmp50, "multi"))
#'     system.time(y <- adiv_table(hmp50, "multi"))
#'     identical(x, y)
#'
with_cache <- function (env, dots, expr) {
  
  args <- lapply(c(as.list(env), dots), eval, env)
  env  <- list2env(args)
  
  
  #________________________________________________________
  # Fetch the maximum cache size.
  #________________________________________________________
  cache_size <- suppressWarnings(abs(as.integer(getOption(
    x       = "rbiom.cache_size", 
    default = Sys.getenv("RBIOM_CACHE_SIZE", unset = "") ))))
  
  if (is.na(cache_size) || !is_scalar_integerish(cache_size))
    cache_size <- 200 * 1024 ^ 2 # 200 MB
  
  if (!isTRUE(cache_size > 0)) return (eval(expr, env))
  
  
  #________________________________________________________
  # Fetch the cache directory.
  #________________________________________________________
  cache <- getOption("rbiom.cache_dir", default = "")
  if (identical(cache, "")) cache <- Sys.getenv("RBIOM_CACHE_DIR", unset = "")
  if (identical(cache, "")) cache <- file.path(tempdir(), "rbiom", "cache")
  if (!is_scalar_character(cache) || !nzchar(cache))              return (eval(expr, env))
  if (!dir.exists(cache) && !dir.create(cache, recursive = TRUE)) return (eval(expr, env))
  
  
  #________________________________________________________
  # Fetch the hashing function.
  #________________________________________________________
  hash <- getOption("rbiom.cache_hash", default = "")
  if (!is.function(hash)) hash <- rlang::hash
  
  
  #________________________________________________________
  # Hash the call into a filename.
  #________________________________________________________
  fn  <- rlang::call_name(rlang::caller_call())
  key <- hash(c(list(fn), args[order(names(args))]))
  fp  <- file.path(cache, paste0(key, ".rds"))
  fp  <- normalizePath(fp, winslash = "/", mustWork = FALSE)
  
  
  #________________________________________________________
  # If the file already exists, just update the timestamp.
  #________________________________________________________
  if (Sys.setFileTime(fp, Sys.time())) {
    
    expr <- readRDS(fp)
    
  } else {
    
    
    #________________________________________________________
    # Don't write objects that are larger than the cache size.
    #________________________________________________________
    expr      <- eval(expr, env)
    expr_size <- as.integer(object.size(expr))
    
    if (expr_size < cache_size) {
      
      #________________________________________________________
      # Delete old files to stay under the maximum.
      #________________________________________________________
      files <- list.files(cache, full.names = TRUE, recursive = TRUE)
      files <- files[rev(order(file.mtime(files)))]
      files <- files[cumsum(file.size(files)) > cache_size - expr_size]
      file.remove(files)
      
      
      saveRDS(expr, file = fp)
    }
    
  }
  
  return (expr)
}
