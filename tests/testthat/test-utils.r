test_that("utils", {
  
  skip_on_cran()
  
  expect_error(eval_envir(new.env(), x = 1, y = 1))
  expect_true(eq(structure(NA, hash = 1), structure(NA, hash = 2)))
  expect_silent(drop_cols(data.frame()))
  expect_null(keep_cols(data.frame()))
  expect_identical(coan(NULL), '')
  expect_identical(vw(NULL), '<none>')
  expect_identical(vw('a'), 'a')
  expect_identical(vw(c('a', 'b')), 'a and b')
  expect_error(require_package('doesnotexist', 'test'))
  expect_error(require_package('doesnotexist', 'test'))
})
