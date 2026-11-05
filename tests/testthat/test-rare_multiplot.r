
# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)))
  env   <- new.env(parent = globalenv())
  r     <- eval(parse(text = code), envir = env)
  built <- function (x) lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
  expect_true(inherits(r, 'patchwork'))
  expect_equal(length(r$patches$plots), length(p$patches$plots))
  for (i in seq_along(p$patches$plots))
    expect_equal(built(r$patches$plots[[i]]), built(p$patches$plots[[i]]), tolerance = 1e-6)
  expect_equal(built(r), built(p), tolerance = 1e-6)
  invisible(r)
}

# A fresh (unmodified) object: other test files edit the shared fixtures in
# place, which adds non-reproducible entries to their history.
fresh5 <- as_rbiom(system.file("extdata", "hmp50.bz2", package = "rbiom"))[1:5]


test_that("rare_multiplot", {

  p <- expect_silent(rare_multiplot(hmp5, stat.by = 'Body Site'))
  expect_s3_class(p, 'patchwork')
  expect_s3_class(p, 'rbiom_plot')

  skip_on_cran()

  expect_identical(as.character(p$cmd), 'rare_multiplot(biom, stat.by = "Body Site")')
  expect_true(grepl("p1 <- ggplot(d1", p$code, fixed = TRUE))
  expect_true(grepl("p2 <- ggplot(d2", p$code, fixed = TRUE))
  expect_true(grepl("patchwork::wrap_plots(p1, p2, ncol = 1)", p$code, fixed = TRUE))
  expect_false(grepl("d1 <- adiv_table", p$code, fixed = TRUE))

  # Both panels use the same rarefaction depth.
  depth <- suggest_rarefy_depth(hmp5$counts)
  expect_true(grepl(sprintf("xintercept = %s", depth), p$code, fixed = TRUE))
  expect_true(grepl(sprintf("yintercept = %s", depth), p$code, fixed = TRUE))

  expect_silent(rare_multiplot(hmp5, rline = FALSE, layers = "p"))
  expect_null(rare_multiplot(hmp5, rline = 300)$stats)
  expect_error(rare_multiplot(hmp5, size = 2), "not an argument")
})


test_that("rare_multiplot code replays", {

  skip_on_cran()

  expect_replays(rare_multiplot(rare50[1:12], stat.by = "Body Site"))
  expect_replays(rare_multiplot(fresh5, rline = 200))
})
