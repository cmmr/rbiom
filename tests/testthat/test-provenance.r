test_that("provenance", {
  
  skip_on_cran()
  
  env <- new.env()
  
  expect_silent(cmd_wrap('doesnotexist', 'dne', env))
  expect_error(env$dne())
  expect_silent(cmd_wrap('base', 'sum', env))
  expect_silent(env$sum(1:4, .lhs = 'x'))
  expect_silent(basewrap('base', 'sum', env))
  expect_silent(env$.sum(1:4))
})
