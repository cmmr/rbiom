test_that("stats_emmeans", {
  
  df <- expect_silent(adiv_table(rare50))
  
  expect_silent(stats_emmeans(df))
  expect_silent(stats_emtrends(df, regr = "Age"))
  
  skip_on_cran()
  
  expect_silent(stats_emmeans(df, regr = "Age"))
  expect_silent(stats_emmeans(df, stat.by = "Sex"))
  expect_silent(stats_emmeans(df, regr = "Age", stat.by = "Sex"))
  expect_silent(stats_emmeans(df, regr = "Age", split.by = "Sex"))
  expect_silent(stats_emmeans(df, regr = "Age", fit = "gam", at = c(30, 90, 150)))
  
  expect_silent(stats_emtrends(df, regr = "Age", stat.by = "Sex"))
  expect_silent(stats_emtrends(df, regr = "Age", split.by = "Sex"))
  expect_silent(stats_emtrends(df, regr = "Age", fit = "gam", at = c(30, 90, 150)))
  
})
