test_that("boxplot", {
  
  df <- expect_silent(adiv_table(rare5))
  
  expect_silent(stats_boxplot(df))
  expect_silent(adiv_boxplot(rare5))
  expect_silent(bdiv_boxplot(rare5))
  expect_silent(taxa_boxplot(rare5))
  
  skip_on_cran()
  
  expect_silent(bdiv_boxplot(rare5, x = 'Body Site', flip = TRUE))
  expect_silent(bdiv_boxplot(rare5, x = 'Body Site', within = 'Body Site'))
  
  expect_silent(taxa_boxplot(rare5, x = "Body Site", rank = c("Phylum", "Genus")))
  expect_silent(taxa_boxplot(rare50, layers = "bed", stat.by = 'Body Site', p.top = 3))
  
  expect_silent(adiv_boxplot(rare50, x = 'Body Site', stat.by = 'Body Site'))
  expect_silent(adiv_boxplot(rare50, x = 'Body Site', stat.by = 'Sex'))
  expect_silent(adiv_boxplot(rare50, stat.by = 'Sex'))
  expect_silent(adiv_boxplot(rare50, facet.by = 'Sex', stat.by = 'Body Site', scales = 'free'))
})
