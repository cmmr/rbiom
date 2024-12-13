test_that("boxplot_build", {
  
  skip_on_cran()
  
  params <- list2env(list(df = structure(list(metadata = NA), class = 'rbiom')))
  expect_error(boxplot_build(params))
  
  df <- adiv_table(rare50)
  expect_warning(stats_boxplot(df, x = 'Sex', facet.by = 'Sex'))
  expect_warning(stats_boxplot(df, stat.by = 'Sex', facet.by = 'Sex'))
  
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'vbx', outliers = FALSE))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', patterns = TRUE, layers = 'dbsl', transform = 'log1p'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'crossbar',   ci = 'range', scales = 'free'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'errorbar',   ci = 'mad',   scales = 'free_x'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'linerange',  ci = 'sd',    scales = 'free_y'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', layers = 'pointrange', ci = 'se',    scales = 'fixed'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', patterns = TRUE, layers = 'cv'))
  expect_silent(stats_boxplot(df, x = 'Sex', stat.by = 'Body Site', flip = TRUE))
})
