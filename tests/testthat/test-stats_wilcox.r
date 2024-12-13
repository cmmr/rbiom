test_that("stats_wilcox", {
  
  skip_on_cran()
  
  df <- expect_silent(adiv_table(rare50))
  expect_silent(stats_wilcox(df))
  expect_error(stats_kruskal(df, stat.by = NULL))
  
})
