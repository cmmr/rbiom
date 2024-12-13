test_that("embed", {
  
  expect_silent(embed_csv(head(mtcars)))
  expect_silent(embed_code('head(mtcars)'))
  
})
