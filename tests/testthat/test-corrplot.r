test_that("corrplot", {
  
  df <- expect_silent(adiv_table(rare5))
  
  expect_silent(stats_corrplot(df,   'age'))
  expect_silent(adiv_corrplot(rare5, 'age'))
  expect_silent(bdiv_corrplot(rare5, 'age'))
  expect_silent(rare_corrplot(hmp5))
  expect_silent(taxa_corrplot(rare5, 'age'))
  
  skip_on_cran()
  
  expect_silent(adiv_corrplot(
    biom      = rare50, 
    x         = 'age', 
    check     = TRUE, 
    stat.by   = 'Sex', 
    layers    = 'ptr', 
    colors    = 'okabe', 
    transform = 'rank' ))
})
