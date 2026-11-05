test_that("cache", {

  skip_on_cran()

  expect_null(rlang::with_options(get_cache_dir(),  rbiom.cache_dir = 'FALSE'))
  expect_null(rlang::with_options(get_cache_dir(),  rbiom.cache_dir = FALSE))
  expect_null(rlang::with_options(get_hash_fun()(), rbiom.cache_dir = 'FALSE'))
  expect_null(rlang::with_options(get_cache_file(), rbiom.cache_dir = 'FALSE'))

  # With caching off, set_cache_value() still returns the result.
  expect_identical(rlang::with_options(set_cache_value(NULL, NA), rbiom.cache_dir = 'FALSE'), NA)
  expect_null(get_cache_value(NULL))

  file.create(f <- tempfile())
  expect_warning(rlang::with_options(get_cache_dir(), rbiom.cache_dir = f))
  unlink(f)

  expect_identical(rlang::hash, rlang::with_options(get_hash_fun(), rbiom.cache_hash = 'FALSE'))
  expect_identical(identity, rlang::with_options(get_hash_fun(), rbiom.cache_hash = identity))
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
