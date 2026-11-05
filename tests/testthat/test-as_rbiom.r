test_that("as_rbiom", {
  old_hmp5 <- as_rbiom(as.list(hmp5), pkg_version = "1")
  expect_equal_rbiom(as_rbiom(hmp5), hmp5)
  expect_equal_rbiom(as_rbiom(old_hmp5), hmp5)
  expect_equal_rbiom(as_rbiom(hmp5$counts), min5)
  expect_equal_rbiom(as_rbiom(as.matrix(hmp5$counts)), min5)
  expect_equal_rbiom(as_rbiom(as.data.frame(as.matrix(hmp5$counts))), min5)
  expect_error(as_rbiom(NULL))
  expect_error(as_rbiom(NA))

  # Same object is returned when nothing changes; a clone otherwise.
  expect_identical(as_rbiom(hmp5), hmp5)
  x <- as_rbiom(hmp5, id = "renamed")
  expect_identical(x$id, "renamed")
  expect_false(identical(x, hmp5))
  expect_error(as_rbiom(hmp5, bogus = 1))
})

test_that("as_rbiom.default", {
  json <- write_biom(hmp5, file = NULL, format = "json")
  expect_equal_rbiom(as_rbiom(json), hmp5)

  df <- as.data.frame(matrix(1:4, 2, dimnames = list(c('a','b'), c('x','y'))))
  expect_s3_class(as_rbiom(df), 'rbiom')
  expect_false(as_rbiom(df)$history[[1]] %>% identical(NULL))
})


test_that("as_rbiom accepts triangular, diagonal, and 1x1 matrices", {
  ids   <- list(c("OTU1", "OTU2", "OTU3"), c("S1", "S2", "S3"))
  cases <- list(
    triangular = matrix(c(600, 300, 600, 0, 300, 500, 0, 0, 7), 3, dimnames = ids),
    diagonal   = matrix(c(5, 0, 0, 0, 8, 0, 0, 0, 4), 3, dimnames = ids),
    one_by_one = matrix(1500, 1, dimnames = list("OTU1", "S1")) )

  for (mtx in cases) {
    biom <- expect_silent(as_rbiom(mtx))
    expect_s4_class(biom$counts, "dgCMatrix")
    expect_equal(as.matrix(biom$counts), mtx)

    path <- tempfile(fileext = ".tsv")
    write_biom(biom, path, format = "tab")
    expect_equal(as.matrix(read_biom(path)$counts), mtx)
    unlink(path)
  }
})
