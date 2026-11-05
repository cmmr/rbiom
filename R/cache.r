#________________________________________________________
# Result caching.
#
# Expensive results are written to `<cache_dir>/<key>.rds`, where `key` is a
# hash of the function name, its (evaluated) arguments, and the versions of
# the packages that determine the result. rbiom objects contribute their
# content hash (`biom$hash`), so two objects with identical contents share
# cache entries regardless of how they were produced. Provenance attributes
# (`cmd`, `history`) are therefore attached *after* a cache lookup - see
# `prov_attach()` in provenance.r.
#________________________________________________________

CACHE <- new.env()


#________________________________________________________
# Fetch the cache directory, or NULL when caching is disabled.
#________________________________________________________
get_cache_dir <- function () {

  cache <- getOption("rbiom.cache_dir", default = "")
  if (isFALSE(cache)) return (NULL)
  if (eq(cache, "")) cache <- Sys.getenv("RBIOM_CACHE_DIR", unset = "")
  if (eq(cache, "")) cache <- file.path(tempdir(), "rbiom", "cache")
  if (!is_scalar_character(cache) || eq(toupper(cache), "FALSE")) return (NULL)
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
# Package versions that are part of every cache key, so that a persistent
# cache is invalidated by upgrades.
#________________________________________________________
cache_versions <- function () {

  if (is.null(CACHE$versions))
    CACHE$versions <- vapply(
      X         = c('rbiom', 'ecodive', 'ggplot2', 'emmeans', 'vegan', 'mgcv'),
      FUN.VALUE = '',
      FUN       = function (p) as.character(packageVersion(p)) )

  return (CACHE$versions)
}


#________________________________________________________
# Hash a function call to a cache file.
# Returns NULL when caching is disabled. Otherwise returns the file path with
# an 'exists' attribute indicating whether a cached result is available.
#________________________________________________________
get_cache_file <- function (fn, params) {

  cache_dir <- get_cache_dir()
  if (is.null(cache_dir)) return (NULL)

  hash <- getOption("rbiom.cache_hash", default = "")
  if (!is.function(hash)) hash <- rlang::hash

  if (is.environment(params)) params <- as.list(params, all.names = TRUE)

  # `cpus` never changes a result.
  params <- params[setdiff(names(params), 'cpus')]
  params <- lapply(params, function (x) { if (inherits(x, "rbiom")) x$hash else hash(x) })

  cache_key  <- hash(list(fn, params[order(names(params))], cache_versions()))
  cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
  cache_file <- normalizePath(cache_file, winslash = "/", mustWork = FALSE)

  # Touch the file so that recently used results are evicted last.
  attr(cache_file, 'exists') <- Sys.setFileTime(cache_file, Sys.time())

  return (cache_file)
}


#________________________________________________________
# Read a cached result, or NULL if there isn't one.
#________________________________________________________
get_cache_value <- function (cache_file) {

  if (!isTRUE(attr(cache_file, 'exists', exact = TRUE))) return (NULL)

  result <- tryCatch(readRDS(cache_file), error = function (e) NULL)

  return (result)
}


#________________________________________________________
# Write a result to the cache. Always returns `result` (visibly), so that
# `set_cache_value(cache_file, x)` can be a function's final expression.
#________________________________________________________
set_cache_value <- function (cache_file, result) {

  if (is.null(cache_file)) return (result)


  #________________________________________________________
  # Fetch the maximum cache size.
  #________________________________________________________
  cache_size <- suppressWarnings(abs(as.numeric(getOption(
    x       = "rbiom.cache_size",
    default = Sys.getenv("RBIOM_CACHE_SIZE", unset = "") ))))

  if (length(cache_size) != 1 || is.na(cache_size))
    cache_size <- 200 * 1024 ^ 2 # 200 MB


  #________________________________________________________
  # Don't write objects that are larger than the cache size.
  #________________________________________________________
  result_size <- as.numeric(object.size(result))

  if (result_size < cache_size) {

    #________________________________________________________
    # Delete the least recently used cache files to stay under
    # the maximum. Only files that rbiom itself wrote are removed.
    #________________________________________________________
    cache <- dirname(cache_file)
    files <- list.files(cache, pattern = '^[0-9a-f]{32}\\.rds$', full.names = TRUE)
    files <- setdiff(files, cache_file)
    files <- files[rev(order(file.mtime(files)))]
    files <- files[cumsum(file.size(files)) > cache_size - result_size]
    if (length(files) > 0) suppressWarnings(file.remove(files))

    # A failed cache write must never fail the user's computation.
    tryCatch(
      expr    = saveRDS(result, file = cache_file),
      error   = function (e) suppressWarnings(file.remove(cache_file)),
      warning = function (w) suppressWarnings(file.remove(cache_file)) )
  }


  return (result)
}
