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
#
# Several R processes may share one cache directory (e.g. a pool of workers
# with the same `rbiom.cache_dir`). So that this is safe and cheap:
#
#   * An entry is written to a temporary file in the cache directory and
#     renamed into place (`cache_write()`), so a reader never finds part of
#     an entry under its final name. The same key always holds the same
#     result, so when two processes write the same entry, the first one
#     to finish is kept.
#   * The directory is listed only now and then, not on every write
#     (`cache_store()`), and least recently used entries are evicted
#     down to 80% of the limit (`cache_evict()`).
#   * Files that other processes delete or replace in the meantime are
#     tolerated everywhere, and a failure to read, write, or evict is a
#     cache miss, never an error or a warning.
#________________________________________________________

CACHE <- new.env()
CACHE$tally <- list() # Per cache directory; see cache_store().

# File names: <key>.rds for entries, <key>.rds.<pid>.<random>.tmp while
# an entry is being written. Only files with these names are ever deleted.
CACHE_ENTRY   <- '^[0-9a-f]{32}\\.rds$'
CACHE_TMP     <- '^[0-9a-f]{32}\\.rds\\.[0-9]+\\.[0-9a-f]{8}\\.tmp$'

CACHE_TARGET  <- 0.8   # Eviction shrinks the cache to this fraction of the limit.
CACHE_RESCAN  <- 100L  # List the directory at least once every this many writes.
CACHE_TMP_AGE <- 3600  # Seconds before an abandoned temporary file is removed.


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
# The maximum size of the cache directory, in bytes.
#________________________________________________________
cache_limit <- function () {

  cache_size <- suppressWarnings(abs(as.numeric(getOption(
    x       = "rbiom.cache_size",
    default = Sys.getenv("RBIOM_CACHE_SIZE", unset = "") ))))

  if (length(cache_size) != 1 || is.na(cache_size))
    cache_size <- 50 * 1024 ^ 2 # 50 MB

  return (cache_size)
}


#________________________________________________________
# Hash a function call to a cache file.
# Returns NULL when caching is disabled. Otherwise returns the file path with
# an 'exists' attribute indicating whether a cached result is available.
#
# Keys are always rlang::hash() (128-bit XXH3, 32 hex digits), which is also
# what `biom$hash` uses, and what CACHE_ENTRY expects when evicting.
#________________________________________________________
get_cache_file <- function (fn, params) {

  cache_dir <- get_cache_dir()
  if (is.null(cache_dir)) return (NULL)

  if (is.environment(params)) params <- as.list(params, all.names = TRUE)

  # `cpus` never changes a result.
  params <- params[setdiff(names(params), 'cpus')]
  params <- lapply(params, function (x) { if (inherits(x, "rbiom")) x$hash else rlang::hash(x) })

  cache_key  <- rlang::hash(list(fn, params[order(names(params))], cache_versions()))
  cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
  cache_file <- normalizePath(cache_file, winslash = "/", mustWork = FALSE)

  return (cache_touch(cache_file))
}


#________________________________________________________
# Set a cache file's modification time to now, so that recently used
# results are evicted last, and record in its 'exists' attribute whether
# there is an entry to read. An entry that can't be touched (e.g. in a
# read-only directory, or briefly locked by another process on Windows)
# can still be read.
#________________________________________________________
cache_touch <- function (cache_file) {

  attr(cache_file, 'exists') <-
    Sys.setFileTime(cache_file, Sys.time()) || file.exists(cache_file)

  return (cache_file)
}


#________________________________________________________
# Read a cached result, or NULL if there isn't one. An entry that is
# unreadable, or that another process evicts between `get_cache_file()` and
# this read, is a miss: the caller recomputes and rewrites it.
#________________________________________________________
get_cache_value <- function (cache_file) {

  if (!isTRUE(attr(cache_file, 'exists', exact = TRUE))) return (NULL)

  read <- function () {
    tryCatch(
      expr  = list(suppressWarnings(readRDS(cache_file))),
      error = function (e) NULL )
  }

  # On Windows, another process touching the entry holds it exclusively
  # for a moment, so try twice.
  result <- read()
  if (is.null(result) && file.exists(cache_file)) {
    Sys.sleep(0.01)
    result <- read()
  }

  # Delete an unreadable entry, so that the caller's set_cache_value()
  # replaces it (cache_write() never overwrites an entry). This silently
  # does nothing if the file is already gone, or is open elsewhere on
  # Windows.
  if (is.null(result)) {
    unlink(cache_file, expand = FALSE)
    return (NULL)
  }

  return (result[[1]])
}


#________________________________________________________
# Write a result to the cache. Always returns `result` (visibly), so that
# `set_cache_value(cache_file, x)` can be a function's final expression.
#________________________________________________________
set_cache_value <- function (cache_file, result) {

  # Callers pass the computation itself, e.g. `set_cache_value(f, build())`.
  # Its errors and warnings are the user's, so evaluate it first, outside
  # the handlers below.
  force(result)

  if (is.null(cache_file)) return (result)

  # A failed write or eviction must never fail the computation.
  tryCatch(
    expr  = suppressWarnings(cache_store(cache_file, result)),
    error = function (e) NULL )

  return (result)
}


#________________________________________________________
# The body of set_cache_value(), which may fail.
#
# Eviction needs the size of the cache directory, but listing and stat-ing
# every entry on every write is O(entries) and would be repeated by every
# process sharing the directory. Instead each process keeps a running
# estimate per directory: the size found at its last listing plus the bytes
# it has written since. It lists the directory (and evicts, if needed) only
# on its first write, when the estimate passes the limit, when it has
# written more than 20% of the limit since its last listing, or after
# CACHE_RESCAN writes.
#
# The estimate doesn't see other processes' writes. But each of P processes
# writes less than 20% of the limit plus one entry between its listings,
# and any listing evicts the whole directory, so the directory can exceed
# the limit by at most about P x (20% of the limit + the largest entry).
#________________________________________________________
cache_store <- function (cache_file, result) {

  limit <- cache_limit()

  # Don't write objects that are larger than the cache size.
  if (as.numeric(utils::object.size(result)) >= limit) return (invisible())

  bytes <- cache_write(cache_file, result)
  if (is.na(bytes)) return (invisible())

  cache_dir <- dirname(cache_file)
  tally     <- CACHE$tally[[cache_dir]]
  if (is.null(tally)) tally <- list(total = NA, since = 0, writes = 0L)

  tally$since  <- tally$since  + bytes
  tally$writes <- tally$writes + 1L

  if ( is.na(tally$total) ||
       tally$total + tally$since > limit ||
       tally$since > (1 - CACHE_TARGET) * limit ||
       tally$writes >= CACHE_RESCAN ) {

    # Reset first, so that if eviction fails it is retried after more
    # writes rather than on every write.
    tally <- list(total = 0, since = 0, writes = 0L)
    CACHE$tally[[cache_dir]] <- tally
    tally$total <- cache_evict(cache_dir, limit, keep = basename(cache_file))
  }

  CACHE$tally[[cache_dir]] <- tally

  return (invisible())
}


#________________________________________________________
# Write `result` to `cache_file` atomically: save it to a temporary file in
# the same directory, then rename that into place. Returns the number of
# bytes added to the cache, or NA if nothing was written. The temporary
# file never outlives this call (unless the process dies).
#
# `save` is a seam for testing.
#________________________________________________________
cache_write <- function (cache_file, result, save = saveRDS) {

  cache_file <- as.character(cache_file)

  # Another process may have written this entry since our cache miss.
  if (file.exists(cache_file)) return (0)

  # Unique per process (pid), and per call and host (random suffix). The
  # suffix doesn't use R's random number generator, so the caller's RNG
  # state is left untouched.
  CACHE$n_tmp <- if (is.null(CACHE$n_tmp)) 1 else CACHE$n_tmp + 1
  suffix      <- substr(rlang::hash(list(
    Sys.info()[['nodename']], Sys.getpid(), Sys.time(), CACHE$n_tmp, tempfile() )), 1, 8)
  tmp <- sprintf("%s.%d.%s.tmp", cache_file, Sys.getpid(), suffix)

  on.exit(unlink(tmp, expand = FALSE))

  saved <- tryCatch(
    expr    = { save(result, tmp); TRUE },
    error   = function (e) FALSE,
    warning = function (w) FALSE ) # e.g. a full disk
  if (!saved) return (NA_real_)

  bytes <- file.size(tmp)
  if (is.na(bytes)) return (NA_real_)


  #________________________________________________________
  # A rename within one directory is atomic: readers find either no entry
  # or a complete one. If another process wrote this entry while we were
  # saving ours, keep theirs - the same key always holds the same result.
  # Besides saving a rename, that matters on Windows, where renaming onto
  # a file that another process has open (e.g. is reading) fails, and
  # file.rename() itself keeps retrying for about 5 seconds first. So a
  # failed rename is retried only briefly, and only while there is still
  # no entry.
  #________________________________________________________
  for (wait in c(0, 0.05)) {
    Sys.sleep(wait)
    if (file.exists(cache_file)) return (0)
    if (suppressWarnings(file.rename(tmp, cache_file))) return (bytes)
    if (!file.exists(tmp)) break
  }

  return (if (file.exists(cache_file)) 0 else NA_real_)
}


#________________________________________________________
# List the cache directory's entries and temporary files, with their sizes
# and modification times. Only names that rbiom itself writes are returned,
# so other files in the directory are never touched. Files deleted by
# another process during the listing are dropped.
#________________________________________________________
cache_scan <- function (cache_dir) {

  name <- list.files(cache_dir)
  tmp  <- grepl(CACHE_TMP, name)
  ok   <- tmp | grepl(CACHE_ENTRY, name)
  name <- name[ok]
  tmp  <- tmp[ok]

  info <- file.info(file.path(cache_dir, name), extra_cols = FALSE)
  ok   <- !is.na(info$size) & !info$isdir

  data.frame(
    name  = name[ok],
    tmp   = tmp[ok],
    size  = info$size[ok],
    mtime = info$mtime[ok] )
}


#________________________________________________________
# Remove temporary files abandoned by processes that died mid-write. Then,
# if the entries total more than `limit` bytes, delete the least recently
# used ones (by modification time, which a cache hit refreshes) until they
# total at most CACHE_TARGET * limit. `keep` names an entry that is never
# deleted - the one just written. Returns the size of the remaining entries.
#
# Another process may be evicting at the same time: deleting a file that
# is already gone, or that is open elsewhere on Windows, fails silently.
#________________________________________________________
cache_evict <- function (cache_dir, limit, keep = NULL) {

  files <- cache_scan(cache_dir)
  del   <- function (x) unlink(file.path(cache_dir, x), expand = FALSE)

  stale <- files$tmp & files$mtime < Sys.time() - CACHE_TMP_AGE
  del(files$name[stale])

  files <- files[!files$tmp,,drop=FALSE]
  total <- sum(files$size)
  if (total <= limit) return (total)

  kept  <- files$name %in% keep
  files <- files[order(!kept, -as.numeric(files$mtime)),,drop=FALSE]
  kept  <- files$name %in% keep
  drop  <- !kept & cumsum(files$size) > CACHE_TARGET * limit
  del(files$name[drop])

  return (total - sum(files$size[drop]))
}
