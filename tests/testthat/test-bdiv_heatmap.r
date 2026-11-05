
# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)))
  env   <- new.env(parent = globalenv())
  r     <- eval(parse(text = code), envir = env)
  built <- function (x) lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
  if (inherits(p, 'patchwork')) {
    expect_true(inherits(r, 'patchwork'))
    expect_equal(length(r$patches$plots), length(p$patches$plots))
    for (i in seq_along(p$patches$plots))
      expect_equal(built(r$patches$plots[[i]]), built(p$patches$plots[[i]]), tolerance = 1e-6)
  }
  expect_equal(built(r), built(p), tolerance = 1e-6)
  invisible(r)
}


test_that("bdiv_heatmap", {

  expect_silent(bdiv_heatmap(rare5, bdiv = c('bray', 'soren'), tracks = c("Body Site", "Age")))

  skip_on_cran()

  expect_silent(bdiv_heatmap(rare5[1:2]))
  expect_silent(bdiv_heatmap(rare5, labs.subtitle = "A subtitle", legend = "bottom"))
  expect_silent(bdiv_heatmap(rare5, tracks = list('Body Site' = list(side = "left")), trees = "cols", label = "rows"))
  expect_error(bdiv_heatmap(rare5[1]), "two samples")
  expect_error(bdiv_heatmap(rare5, title = NA))
  expect_error(bdiv_heatmap(rare5, subtitle = "x"), "labs.")

  p <- bdiv_heatmap(rare5, bdiv = "jaccard")
  expect_identical(as.character(p$cmd), 'bdiv_heatmap(biom, "jaccard")')
  expect_true(grepl("Jaccard", p$code, fixed = TRUE))
  expect_true(grepl("dm  <- bdiv_distmat(biom", p$code, fixed = TRUE))

  # Samples are clustered on the plotted distances.
  hc <- stats::hclust(bdiv_distmat(rare5, "jaccard"), method = "complete")
  expect_identical(levels(p$data$x), hc$labels[hc$order])
})


test_that("bdiv_heatmap code replays", {

  skip_on_cran()

  b <- rare50[1:12]
  expect_replays(bdiv_heatmap(b, bdiv = "bray", tracks = "Sex"))
  expect_replays(bdiv_heatmap(b, bdiv = c("bray", "jaccard"), tracks = "Body Site"))
  expect_replays(bdiv_heatmap(b, tracks = list('Sex' = list(side = "left"), 'Age' = "reds"), trees = FALSE))
  expect_replays(bdiv_heatmap(rare5[1:2]))
})
