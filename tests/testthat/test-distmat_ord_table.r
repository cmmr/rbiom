test_that("distmat_ord_table", {
  
  dm15 <- expect_silent(bdiv_distmat(rare50[1:15]))
  expect_silent(distmat_ord_table(dm15, ord = "PCoA"))
  
  skip_on_cran()
  
  dm40 <- expect_silent(bdiv_distmat(rare50[1:40]))
  expect_silent(distmat_ord_table(dm15, ord = "PCoA"))
  expect_silent(distmat_ord_table(dm15, ord = "tSNE"))
  expect_silent(distmat_ord_table(dm15, ord = "UMAP"))
  expect_silent(distmat_ord_table(dm40, ord = "NMDS"))
  
})
