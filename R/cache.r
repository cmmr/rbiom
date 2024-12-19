

#________________________________________________________
# Fetch the cache directory.
#________________________________________________________
get_cache_dir <- function () {
  
  cache <- getOption("rbiom.cache_dir", default = "")
  if (eq(cache, "")) cache <- Sys.getenv("RBIOM_CACHE_DIR", unset = "")
  if (eq(cache, "")) cache <- file.path(tempdir(), "rbiom", "cache")
  if (!is_scalar_character(cache) || eq(cache, "FALSE")) return (NULL)
  if (!dir.exists(cache) && !dir.create(cache, recursive = TRUE)) return (NULL)
  
  return (normalizePath(cache, winslash = '/'))
}


#________________________________________________________
# The function for hashing cached objects.
#________________________________________________________
get_hash_fun <- function () {
  
  if (is.null(get_cache_dir())) {
    hash_fun <- function (obj) NULL
    
  } else {
    hash_fun <- getOption("rbiom.cache_hash", default = "")
    if (!is.function(hash_fun)) hash_fun <- rlang::hash
  }
  
  return (hash_fun)
}



#________________________________________________________
# Hash a function call to a cache file.
#________________________________________________________
get_cache_file <- function (fn, params) {
  
  cache_dir <- get_cache_dir()
  if (is.null(cache_dir)) return (NULL)
  
  hash <- getOption("rbiom.cache_hash", default = "")
  if (!is.function(hash)) hash <- rlang::hash
  
  params  %<>% lapply(function (x) { if (inherits(x, "rbiom")) x$hash else hash(x) })
  cache_key  <- hash(c(list(fn), params[order(names(params))]))
  cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
  cache_file <- normalizePath(cache_file, winslash = "/", mustWork = FALSE)
  
  attr(cache_file, 'exists') <- Sys.setFileTime(cache_file, Sys.time())
  
  # cat(sprintf("%s => %s [%s]\n", fn, cache_key, attr(cache_file, 'exists')))
  
  return (cache_file)
}



set_cache_value <- function (cache_file, result) {
  
  if (is.null(cache_file)) return (NULL)
  
  
  #________________________________________________________
  # Fetch the maximum cache size.
  #________________________________________________________
  cache_size <- suppressWarnings(abs(as.integer(getOption(
    x       = "rbiom.cache_size", 
    default = Sys.getenv("RBIOM_CACHE_SIZE", unset = "") ))))
  
  if (is.na(cache_size) || !is_scalar_integerish(cache_size))
    cache_size <- 200 * 1024 ^ 2 # 200 MB
  
  
  
  #________________________________________________________
  # Don't write objects that are larger than the cache size.
  #________________________________________________________
  result_size <- as.integer(object.size(result))
  
  if (result_size < cache_size) {
    
    #________________________________________________________
    # Delete old files to stay under the maximum.
    #________________________________________________________
    cache <- dirname(cache_file)
    files <- list.files(cache, full.names = TRUE, recursive = TRUE)
    files <- files[rev(order(file.mtime(files)))]
    files <- files[cumsum(file.size(files)) > cache_size - result_size]
    file.remove(files)
    
    
    saveRDS(result, file = cache_file)
  }
  
  
  return (NULL)
}
