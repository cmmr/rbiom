
# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)))
  env   <- new.env(parent = globalenv())
  r     <- eval(parse(text = code), envir = env)
  built <- function (x) lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
  expect_equal(built(r), built(p), tolerance = 1e-6)
  invisible(r)
}

# A fresh (unmodified) object: other test files edit the shared fixtures in
# place, which adds non-reproducible entries to their history.
fresh5 <- as_rbiom(system.file("extdata", "hmp50.bz2", package = "rbiom"))[1:5]


test_that("rare_stacked", {

  expect_silent(rare_stacked(hmp5, rline = FALSE, y.transform = 'none'))

  skip_on_cran()

  expect_silent(rare_stacked(hmp5))
  expect_silent(rare_stacked(hmp5, rline = 500, r.linewidth = 2, r.linetype = "twodash"))
  expect_silent(rare_stacked(hmp5, counts = FALSE, labels = FALSE))
  expect_error(rare_stacked(hmp5, rline = NA))
  expect_error(rare_stacked(hmp5, y.transform = "sqrt"))
  expect_error(rare_stacked(hmp5, size = 2), "not an argument")

  # `labels = TRUE` really shows the sample names.
  p <- rare_stacked(hmp5, labels = TRUE)
  expect_true(grepl("geom_point", p$code, fixed = TRUE))
  expect_true(grepl("x += NULL", p$code))
  p <- rare_stacked(hmp5, labels = FALSE)
  expect_false(grepl("geom_point", p$code, fixed = TRUE))

  p <- rare_stacked(hmp5, rline = 500)
  expect_identical(as.character(p$cmd), 'rare_stacked(biom, 500)')
  expect_true(grepl("yintercept = 500", p$code, fixed = TRUE))
  expect_true(grepl("Samples Retained", p$code, fixed = TRUE))
  expect_identical(levels(p$data$.group), c("Excluded", "Retained"))

  p <- rare_stacked(hmp5, rline = FALSE)
  expect_false(hasName(p$data, '.group'))
  expect_false(grepl("geom_hline", p$code, fixed = TRUE))
})


test_that("rare_stacked code replays", {

  skip_on_cran()

  expect_replays(rare_stacked(as_rbiom(system.file("extdata", "hmp50.bz2", package = "rbiom"))))
  expect_replays(rare_stacked(fresh5, rline = FALSE, y.transform = 'none', labels = FALSE))
  expect_replays(rare_stacked(fresh5, rline = 1000, counts = FALSE))
})
