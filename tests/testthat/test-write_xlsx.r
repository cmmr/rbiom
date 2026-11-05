test_that("write_xlsx", {

  skip_on_cran()
  skip_if_not_installed('openxlsx')

  tfile <- tempfile(fileext = '.xlsx')
  on.exit(unlink(tfile), add = TRUE)

  biom <- hmp5$clone()
  attr(biom, 'mtcars') <- head(mtcars)
  attr(biom, 'mtx')    <- as.matrix(head(mtcars))
  attr(biom, 'dm')     <- dist(head(mtcars))

  expect_silent(write_xlsx(biom, file = tfile))
  expect_error(write_xlsx(biom, file = tfile), "overwrite")

  sheets <- openxlsx::getSheetNames(tfile)
  expect_true(all(c('Reads Per Sample', 'Sample Metadata', 'Alpha Diversity', 'OTU Taxonomy', 'mtcars', 'mtx', 'dm', 'Genus') %in% sheets))

  adiv <- openxlsx::read.xlsx(tfile, sheet = 'Alpha Diversity')
  expect_identical(names(adiv)[[1]], '.sample')
  unlink(tfile)

  expect_silent(write_xlsx(rare5, file = tfile))
  unlink(tfile)

  # Non-rbiom input, and relative abundances (no rarefaction sheets).
  expect_silent(write_xlsx(as.matrix(hmp5$counts), file = tfile))
  unlink(tfile)
  expect_silent(write_xlsx(biom_relativize(hmp5), file = tfile, overwrite = TRUE))

})
