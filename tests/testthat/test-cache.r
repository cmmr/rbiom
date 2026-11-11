test_that("cache", {

  skip_on_cran()

  expect_null(rlang::with_options(get_cache_dir(),  rbiom.cache_dir = 'FALSE'))
  expect_null(rlang::with_options(get_cache_dir(),  rbiom.cache_dir = FALSE))
  expect_null(rlang::with_options(get_cache_file(), rbiom.cache_dir = 'FALSE'))

  # With caching off, set_cache_value() still returns the result.
  expect_identical(rlang::with_options(set_cache_value(NULL, NA), rbiom.cache_dir = 'FALSE'), NA)
  expect_null(get_cache_value(NULL))

  file.create(f <- tempfile())
  expect_warning(rlang::with_options(get_cache_dir(), rbiom.cache_dir = f))
  unlink(f)

  # Keys are always rlang::hash(): `rbiom.cache_hash` is no longer used, so
  # every entry has a name that eviction recognizes.
  dir <- tempfile()
  f1  <- rlang::with_options(get_cache_file('x', list(a = 1)), rbiom.cache_dir = dir)
  f2  <- rlang::with_options(get_cache_file('x', list(a = 1)), rbiom.cache_dir = dir, rbiom.cache_hash = identity)
  expect_identical(as.character(f1), as.character(f2))
  expect_match(basename(f1), CACHE_ENTRY)
})


test_that("results are identical with caching on and off", {

  skip_on_cran()

  on  <- rlang::with_options(list(
    a = adiv_table(hmp5),
    b = bdiv_distmat(hmp5, 'bray'),
    t = taxa_table(hmp5, rank = 'Phylum') ), rbiom.cache_dir = tempfile())

  off <- rlang::with_options(list(
    a = adiv_table(hmp5),
    b = bdiv_distmat(hmp5, 'bray'),
    t = taxa_table(hmp5, rank = 'Phylum') ), rbiom.cache_dir = 'FALSE')

  expect_equal(on, off)
  expect_s3_class(off$a, 'tbl_df')
  expect_s3_class(off$b, 'dist')
})


test_that("cache hits and safe eviction", {

  skip_on_cran()

  dir <- tempfile()

  rlang::with_options(rbiom.cache_dir = dir, {

    # A second call is served from the cache.
    f1 <- get_cache_file('adiv_table', list(biom = hmp5))
    expect_false(attr(f1, 'exists'))
    x  <- adiv_table(hmp5)
    f2 <- get_cache_file('adiv_table', list(biom = hmp5))
    expect_identical(as.character(f1), as.character(f2))

    # Only rbiom's own files are ever deleted when trimming the cache.
    user_file <- file.path(get_cache_dir(), 'my_results.rds')
    saveRDS(1:10, user_file)

    rlang::with_options(rbiom.cache_size = 1, {
      set_cache_value(get_cache_file('x', list(a = 1)), 1)
    })
    expect_true(file.exists(user_file))
  })
})



#________________________________________________________
# A fresh cache directory for one test.
#________________________________________________________
local_cache_dir <- function (size = NULL, env = parent.frame()) {
  dir <- tempfile('cache_')
  dir.create(dir)
  dir <- normalizePath(dir, winslash = '/')
  rlang::local_options(rbiom.cache_dir = dir, rbiom.cache_size = size, .frame = env)
  return (dir)
}

# Write an entry directly, as another process would, last used `mins` ago.
put_entry <- function (dir, i, value, mins = 0) {
  path <- file.path(dir, sprintf("%032x.rds", i))
  saveRDS(value, path, compress = FALSE)
  Sys.setFileTime(path, Sys.time() - 60 * mins)
  return (path)
}


test_that("cache size: 50 MB by default; the option overrides the env var", {

  old <- Sys.getenv('RBIOM_CACHE_SIZE', unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv('RBIOM_CACHE_SIZE') else Sys.setenv(RBIOM_CACHE_SIZE = old))

  Sys.unsetenv('RBIOM_CACHE_SIZE')
  rlang::local_options(rbiom.cache_size = NULL)
  expect_equal(cache_limit(), 50 * 1024 ^ 2)

  Sys.setenv(RBIOM_CACHE_SIZE = 1024 ^ 3)
  expect_equal(cache_limit(), 1024 ^ 3)

  rlang::local_options(rbiom.cache_size = 2e6)
  expect_equal(cache_limit(), 2e6)
})


test_that("a cached value is read back; an unreadable entry is a miss", {

  skip_on_cran()
  dir <- local_cache_dir()

  f <- get_cache_file('x', list(a = 1))
  expect_false(attr(f, 'exists'))
  expect_null(get_cache_value(f))
  expect_identical(set_cache_value(f, mtcars), mtcars)

  f <- get_cache_file('x', list(a = 1))
  expect_true(attr(f, 'exists'))
  expect_identical(get_cache_value(f), mtcars)

  # Garbage under an entry's name is a silent miss, and is deleted so that
  # the recomputed result replaces it.
  writeBin(as.raw(1:100), f)
  f <- get_cache_file('x', list(a = 1))
  expect_true(attr(f, 'exists'))
  expect_silent(expect_null(get_cache_value(f)))
  expect_false(file.exists(f))
  set_cache_value(f, mtcars)
  expect_identical(get_cache_value(get_cache_file('x', list(a = 1))), mtcars)

  # So is an entry deleted between lookup and read (e.g. evicted by another
  # process).
  f <- get_cache_file('x', list(a = 1))
  unlink(f)
  expect_silent(expect_null(get_cache_value(f)))

  # Results larger than the cache size are not cached.
  rlang::local_options(rbiom.cache_size = 1000)
  f <- get_cache_file('x', list(a = 2))
  expect_identical(set_cache_value(f, 1:1000 + 0), 1:1000 + 0)
  expect_false(file.exists(f))
})


test_that("an entry appears under its final name only when complete", {

  skip_on_cran()
  dir <- local_cache_dir()
  f   <- get_cache_file('x', list(a = 1))

  # While the value is being written, it is only in a temporary file in
  # the cache directory.
  seen  <- NULL
  bytes <- cache_write(f, mtcars, save = function (x, file) {
    saveRDS(x, file)
    seen <<- list(final = file.exists(f), files = list.files(dir), tmp = basename(file))
  })
  expect_false(seen$final)
  expect_identical(seen$files, seen$tmp)
  expect_match(seen$tmp, CACHE_TMP)
  expect_equal(bytes, file.size(f))
  expect_identical(list.files(dir), basename(f))
  expect_identical(readRDS(f), mtcars)

  # An existing entry is kept (the same key holds the same result).
  expect_identical(cache_write(f, mtcars, save = function (...) stop("unused")), 0)
  expect_identical(readRDS(f), mtcars)

  # A failed write leaves neither an entry nor a temporary file.
  g <- get_cache_file('x', list(a = 2))
  for (fail in c(stop, warning)) {
    expect_identical(NA_real_, cache_write(g, mtcars, save = function (x, file) {
      writeBin(as.raw(1:100), file)
      fail("disk full")
    }))
    expect_identical(list.files(dir), basename(f))
  }

  # The computation's own errors and warnings reach the caller.
  expect_error(set_cache_value(g, stop("bad input")), "bad input")
  expect_warning(x <- set_cache_value(g, { warning("odd input"); 1 }), "odd input")
  expect_identical(x, 1)

  # A failed write doesn't fail the computation.
  local_mocked_bindings(cache_write = function (...) stop("disk full"))
  expect_silent(x <- set_cache_value(g, mtcars))
  expect_identical(x, mtcars)
  local_mocked_bindings(cache_write = function (...) warning("disk full"))
  expect_silent(x <- set_cache_value(g, mtcars))
  expect_identical(x, mtcars)
})


test_that("the cache directory is listed only when needed", {

  skip_on_cran()

  limit <- 1e6
  dir   <- local_cache_dir(size = limit)

  n    <- 0
  scan <- cache_scan
  local_mocked_bindings(cache_scan = function (...) { n <<- n + 1; scan(...) })

  put <- function (i, x) set_cache_value(get_cache_file('x', list(i = i)), x)
  big <- function (k) sin(seq_len(k / 8)) # ~k bytes; doesn't compress well

  # 900 KB written by another process, used an hour ago.
  other <- put_entry(dir, 1, raw(9e5), mins = 60)

  # This process's first write to the directory lists it.
  put(0, 1)
  expect_identical(n, 1)

  # Writes that keep the estimate under the limit don't.
  for (i in 1:20) put(i, i)
  expect_identical(n, 1)

  # A write that passes the limit lists the directory and evicts the least
  # recently used entries, but not the one just written.
  put(21, big(1.5e5))
  expect_identical(n, 2)
  expect_false(file.exists(other))
  expect_true(attr(get_cache_file('x', list(i = 21)), 'exists'))

  # So does writing 20% of the limit since the last listing, which bounds
  # how far several processes can overshoot the limit together.
  put(22, big(1.3e5))
  expect_identical(n, 2)
  put(23, big(1.3e5))
  expect_identical(n, 3)

  # And so do CACHE_RESCAN writes, however small.
  for (i in seq_len(CACHE_RESCAN - 1)) put(100 + i, i)
  expect_identical(n, 3)
  put(99, 1)
  expect_identical(n, 4)

  # A failed listing doesn't fail the computation.
  local_mocked_bindings(cache_scan = function (...) stop("gone"))
  expect_silent(x <- put(24, big(5e5)))
  expect_identical(x, big(5e5))
})


test_that("eviction removes the least recently used entries", {

  skip_on_cran()
  dir <- local_cache_dir()

  # Six ~100 KB entries; the first was used least recently.
  paths <- vapply(1:6, function (i) put_entry(dir, i, as.raw(rep(i, 1e5)), mins = 7 - i), '')

  # Reading the first one makes it the most recently used.
  expect_identical(get_cache_value(cache_touch(paths[[1]])), as.raw(rep(1, 1e5)))

  # Only rbiom's own files are ever deleted.
  user_file <- file.path(dir, 'my_results.rds')
  saveRDS(1:10, user_file)
  Sys.setFileTime(user_file, Sys.time() - 3600)

  # Over the 520 KB limit, so evict down to 80% of it (416 KB).
  size <- cache_evict(dir, limit = 5.2e5)
  expect_identical(file.exists(paths), c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(size, sum(file.size(paths[c(1, 4:6)])))
  expect_true(file.exists(user_file))

  # Under the limit: nothing is deleted.
  expect_equal(cache_evict(dir, limit = 5.2e5), size)
  expect_identical(file.exists(paths), c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE))

  # Files that vanish (e.g. evicted by another process) between listing
  # and deletion are fine.
  scan <- cache_scan
  local_mocked_bindings(cache_scan = function (dir) {
    x <- scan(dir)
    unlink(file.path(dir, x$name))
    x
  })
  expect_silent(cache_evict(dir, limit = 1))
  expect_identical(list.files(dir), basename(user_file))
})


test_that("abandoned temporary files are removed once stale", {

  skip_on_cran()
  dir <- local_cache_dir()

  tmp_file <- function (i, mins) {
    path <- file.path(dir, sprintf("%032x.rds.%d.%08x.tmp", i, 1234L, i))
    writeBin(raw(10), path)
    Sys.setFileTime(path, Sys.time() - 60 * mins)
    return (path)
  }

  stale <- tmp_file(1, mins = 120)
  fresh <- tmp_file(2, mins = 1)
  other <- file.path(dir, 'notes.tmp')
  writeBin(raw(10), other)
  Sys.setFileTime(other, Sys.time() - 7200)

  expect_match(basename(c(stale, fresh)), CACHE_TMP)
  expect_identical(cache_evict(dir, limit = 1e9), 0)
  expect_false(file.exists(stale))
  expect_true(file.exists(fresh))
  expect_true(file.exists(other))
})


test_that("processes sharing a cache directory don't interfere", {

  skip_on_cran()
  skip_if_not_installed("callr")

  dir <- tempfile('cache_')
  dir.create(dir)
  dir <- normalizePath(dir, winslash = '/')
  go  <- file.path(dir, 'go') # Not a cache file name, so never evicted.


  #________________________________________________________
  # The cache functions, detached from the rbiom namespace so that workers
  # run this version of them without loading rbiom (during devtools::test()
  # the installed rbiom may be older).
  #________________________________________________________
  ns  <- asNamespace('rbiom')
  env <- new.env(parent = globalenv())
  for (nm in c(ls(ns, pattern = '^(cache_|CACHE)'), 'get_cache_value', 'set_cache_value')) {
    obj <- get(nm, envir = ns)
    if (is.function(obj)) environment(obj) <- env
    assign(nm, obj, envir = env)
  }
  env$CACHE <- new.env()
  env$CACHE$tally <- list()


  #________________________________________________________
  # Each worker repeatedly looks up 8 keys, computing and caching the ones
  # it doesn't find. Only ~5 entries fit, so entries are evicted and
  # rewritten all the time.
  #________________________________________________________
  worker <- function (env, dir, go, keys, n) {

    options(rbiom.cache_size = 4e5)
    value <- function (k) sin(seq_len(1e4) * k) # 80 KB

    t0 <- Sys.time()
    while (!file.exists(go) && Sys.time() - t0 < 60) Sys.sleep(0.01)

    errors  <- character(0)
    warns   <- character(0)
    corrupt <- character(0)
    wrong   <- 0L
    hits    <- 0L

    for (i in seq_len(n)) {
      k <- keys[[i %% length(keys) + 1]]
      f <- file.path(dir, sprintf("%032x.rds", k))

      # Any entry under its final name must be complete. It may vanish
      # (evicted) or, on Windows, be locked for a moment (touched), so
      # "cannot open" is fine.
      if (file.exists(f)) {
        e <- tryCatch(
          expr  = { suppressWarnings(readRDS(f)); NULL },
          error = function (e) conditionMessage(e) )
        if (!is.null(e) && !grepl('cannot open', e))
          corrupt <- c(corrupt, e)
      }

      withCallingHandlers(
        warning = function (w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart('muffleWarning')
        },
        tryCatch(
          error = function (e) errors <<- c(errors, conditionMessage(e)),
          expr  = {
            x <- env$get_cache_value(env$cache_touch(f))
            if (is.null(x)) {
              env$set_cache_value(f, value(k))
            } else {
              hits <- hits + 1L
              if (!identical(x, value(k))) wrong <- wrong + 1L
            }
          }))
    }

    list(errors = errors, warnings = warns, corrupt = corrupt, wrong = wrong, hits = hits)
  }

  procs <- lapply(1:4, function (w) {
    callr::r_bg(
      func = worker,
      args = list(env = env, dir = dir, go = go, keys = c(w:8, seq_len(w - 1)), n = 200) )
  })
  file.create(go)

  for (p in procs) p$wait(timeout = 120 * 1000)
  for (p in procs) if (p$is_alive()) p$kill()
  res <- lapply(procs, function (p) p$get_result())

  for (r in res) {
    expect_identical(r$errors,   character(0))
    expect_identical(r$warnings, character(0))
    expect_identical(r$corrupt,  character(0))
    expect_identical(r$wrong,    0L)
  }
  expect_gt(sum(vapply(res, `[[`, 0L, 'hits')), 0)

  # Every entry left is complete, and no temporary files are left behind.
  files <- setdiff(list.files(dir), basename(go))
  expect_true(length(files) > 0)
  expect_true(all(grepl(CACHE_ENTRY, files)))
  for (f in files) {
    k <- strtoi(substr(f, 25, 32), 16L)
    expect_identical(readRDS(file.path(dir, f)), sin(seq_len(1e4) * k))
  }
})


test_that("a cache directory that can't be written to doesn't break a plot", {

  skip_on_cran()

  dir <- tempfile('cache_')
  dir.create(dir)

  # Make `dir` read-only.
  if (.Platform$OS.type == 'windows') {
    acl <- function (...) suppressWarnings(system2('icacls', c(shQuote(dir), ...), stdout = FALSE, stderr = FALSE))
    on.exit(acl('/remove:d', '*S-1-1-0'), add = TRUE)
    acl('/deny', '*S-1-1-0:(OI)(CI)(W)') # Everyone
  } else {
    on.exit(Sys.chmod(dir, '0755'), add = TRUE)
    Sys.chmod(dir, '0555')
  }
  if (suppressWarnings(file.create(file.path(dir, 'probe'))))
    skip("Can't make a read-only directory here (e.g. running as root).")

  rlang::local_options(rbiom.cache_dir = dir)
  expect_no_warning(expect_no_error(p <- adiv_boxplot(hmp5)))
  expect_s3_class(p, 'ggplot')
  expect_length(list.files(dir), 0)

  # A cache directory that can't be created only warns, as before.
  file.create(f <- tempfile())
  rlang::local_options(rbiom.cache_dir = file.path(f, 'cache'))
  expect_no_error(p <- suppressWarnings(adiv_boxplot(hmp5)))
  expect_s3_class(p, 'ggplot')
})


test_that("arguments that can't change a result aren't part of its key", {

  skip_on_cran()
  dir <- local_cache_dir()
  n_entries <- function () length(list.files(dir))

  # `underscores` in tree_subset().
  tree <- hmp50$tree
  t1   <- tree_subset(tree, 1:5, underscores = FALSE)
  n    <- n_entries()
  t2   <- tree_subset(tree, 1:5, underscores = TRUE)
  expect_identical(n_entries(), n)
  expect_identical(t1, t2)

  # Provenance of `dm` in distmat_ord_table() and distmat_stats(), except
  # whether it has any: that changes the code in the result.
  dm1 <- bdiv_distmat(hmp50, 'bray')
  dm2 <- dm1
  attr(dm2, 'cmd')     <- 'dm <- other_pipeline()'
  attr(dm2, 'history') <- NULL
  grp <- pull(hmp50, 'Sex')

  o1 <- distmat_ord_table(dm1)
  s1 <- distmat_stats(dm1, grp, permutations = 99)
  n  <- n_entries()
  o2 <- distmat_ord_table(dm2)
  s2 <- distmat_stats(dm2, grp, permutations = 99)
  expect_identical(n_entries(), n)
  expect_equal(strip_prov(o1), strip_prov(o2))
  expect_identical(attr(o2, 'data_cmd'), 'dm <- other_pipeline()')
  expect_identical(attr(s2, 'data_cmd'), 'dm <- other_pipeline()')

  s3 <- distmat_stats(strip_prov(dm1), grp, permutations = 99)
  expect_identical(n_entries(), n + 1L)
  expect_true(any(grepl("not reproducible", s3$code, fixed = TRUE)))
  expect_false(any(grepl("not reproducible", s1$code, fixed = TRUE)))
})
