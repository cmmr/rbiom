test_that("distmat_ord_table", {
  
  dm15 <- bdiv_distmat(rare50[1:15])
  expect_no_error(distmat_ord_table(dm15, ord = "PCoA"))
  
  skip_on_cran()
  
  dm40 <- bdiv_distmat(rare50[1:40])
  expect_no_error(distmat_ord_table(dm15, ord = "PCoA"))
  expect_no_error(distmat_ord_table(dm15, ord = "tSNE"))
  expect_no_error(distmat_ord_table(dm15, ord = "UMAP"))
  expect_no_error(distmat_ord_table(dm40, ord = "NMDS"))
  
})
