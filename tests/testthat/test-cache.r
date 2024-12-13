test_that("cache", {
  
  skip_on_cran()
  
  expect_null(rlang::with_options(get_cache_dir(),  rbiom.cache_dir = 'FALSE'))
  expect_null(rlang::with_options(get_hash_fun()(), rbiom.cache_dir = 'FALSE'))
  expect_null(rlang::with_options(get_cache_file(), rbiom.cache_dir = 'FALSE'))
  expect_null(rlang::with_options(set_cache_value(NULL, NA), rbiom.cache_dir = 'FALSE'))
  
  file.create(f <- tempfile())
  expect_warning(rlang::with_options(get_cache_dir(), rbiom.cache_dir = f))
  unlink(f)
  
  expect_identical(rlang::hash, rlang::with_options(get_hash_fun(), rbiom.cache_hash = 'FALSE'))
  expect_identical(identity, rlang::with_options(get_hash_fun(), rbiom.cache_hash = identity))
  
})
