test_that("bdiv_ord_plot", {
  
  expect_silent(bdiv_ord_plot(rare50))
  expect_silent(bdiv_ord_plot(rare50, stat.by='Body Site'))
  
  skip_on_cran()
  
  expect_silent(bdiv_ord_plot(rare50, layers='sa', stat.by='Body Site', facet.by='Sex'))
  expect_silent(bdiv_ord_plot(rare50, layers='tm', taxa=c('Bacteroides', 'Prevotella')))
  expect_silent(bdiv_ord_plot(rare50, layers='mn', bdiv='Manhattan', weighted=FALSE))
  expect_silent(bdiv_ord_plot(rare50, layers='pt', bdiv='Jaccard',   weighted=FALSE))
  
  expect_silent(bdiv_ord_plot(rare50, rank=NULL))
  expect_silent(bdiv_ord_plot(rare50, rank=c('Phylum', 'Genus')))
  expect_silent(bdiv_ord_plot(rare50, bdiv=c('Manhattan', 'Jaccard')))
  expect_silent(bdiv_ord_plot(rare50, bdiv=c('Manhattan', 'Jaccard'), facet.by='Sex'))
  expect_silent(bdiv_ord_plot(rare50, bdiv=c('Manhattan', 'Jaccard'), rank=c('Phylum', 'Genus')))
  expect_silent(bdiv_ord_plot(rare50, bdiv=c('Manhattan', 'Jaccard'), rank=c('Phylum', 'Genus'), facet.by='Sex'))
  
  expect_error(bdiv_ord_plot(rare50[1:3]))
  
})
