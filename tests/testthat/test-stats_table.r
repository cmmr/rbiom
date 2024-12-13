test_that("stats_table", {
  
  expect_silent(adiv_stats(min5, fit = 'lm'))
  expect_silent(bdiv_stats(min5, fit = 'lm'))
  expect_silent(taxa_stats(hmp5, fit = 'lm', rank = 'p', taxa = 2))
  
})
