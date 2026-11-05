
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

mtx           <- matrix(1:25, ncol = 5)
dimnames(mtx) <- list(letters[1:5], LETTERS[1:5])
cat_vals      <- c("Male", "Female")[c(1,2,1,2,1)]
num_vals      <- c(25, 30, 30, 35, 40)
tracks        <- list(
  list(label = "Sex", values = cat_vals, side = "left", colors   = "bright"),
  list(label = "Age", values = num_vals, side = "left", range    = c(20,40), bins = 3),
  list(label = "Sex", values = cat_vals, side = "top",  na.color = "gray"),
  list(label = "Age", values = num_vals, side = "top",  colors   = c('red', 'blue')) )


test_that("plot_heatmap", {

  expect_silent(plot_heatmap(mtx = mtx, tracks = tracks))

  skip_on_cran()

  expect_silent(plot_heatmap(mtx = mtx, grid = NULL, rescale = 'rows', clust = FALSE))
  expect_silent(plot_heatmap(mtx = mtx, grid = 'greens', rescale = 'cols'))
  expect_silent(plot_heatmap(mtx = mtx, grid = c('green', 'blue')))
  expect_silent(plot_heatmap(mtx = mtx, tracks = list(list(label = "Sex", values = cat_vals, side = "left"))))
  expect_silent(plot_heatmap(mtx = mtx, tracks = list('Sex' = setNames(cat_vals, LETTERS[1:5]))))
  expect_silent(plot_heatmap(mtx = mtx, title = "Title", labs.subtitle = "Sub", legend = "bottom", xlab.angle = 90))
  expect_silent(plot_heatmap(mtx = mtx, dist = stats::dist(mtx), clust = stats::hclust(stats::dist(mtx))))
  expect_silent(plot_heatmap(mtx = mtx[1, , drop = FALSE]))
  expect_silent(plot_heatmap(mtx = mtx[, 1, drop = FALSE]))

  expect_error(plot_heatmap(mtx = unname(mtx)), "names")
  expect_error(plot_heatmap(mtx = mtx, tracks = list(list(label = "x", values = 1:3))), "Length mismatch")
  expect_error(plot_heatmap(mtx = mtx, tracks = list(list(label = "x", values = setNames(1:5, letters[2:6]), side = "left"))), "missing")
  expect_error(plot_heatmap(mtx = mtx, subtitle = "x"), "labs.")

  p <- plot_heatmap(mtx = mtx, grid = "greens")
  expect_identical(as.character(p$cmd), 'plot_heatmap(mtx, "greens")')

  # No vectorised element_text() in the code (it warns on replay).
  p <- plot_heatmap(mtx = mtx, tracks = tracks)
  expect_silent(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p)))
  expect_true(grepl("sec_axis", p$code, fixed = TRUE))
})


test_that("plot_heatmap geometry", {

  skip_on_cran()

  # The second element of tree_height sets the left (row) dendrogram.
  p1 <- plot_heatmap(mtx, tree_height = c(10, 10))
  p2 <- plot_heatmap(mtx, tree_height = c(10, 30))
  t1 <- attr(p1$data, 'trees') %||% ggplot2::layer_data(p1, 3)
  left_x <- function (p) {
    d <- ggplot2::ggplot_build(p)$data
    segs <- d[[which(sapply(d, function (x) all(c("x", "xend") %in% names(x))))[[1]]]]
    min(segs$x, segs$xend)
  }
  expect_true(left_x(p2) < left_x(p1))

  # Left tracks push the row dendrogram further left by the track width.
  p3 <- plot_heatmap(mtx, tracks = tracks[1:2], track_height = c(10, 10))
  p4 <- plot_heatmap(mtx, tracks = tracks[1:2], track_height = c(10, 40))
  expect_true(left_x(p4) < left_x(p3))
  expect_true(left_x(p3) < left_x(p1))
})


test_that("plot_heatmap code replays", {

  skip_on_cran()

  expect_replays(plot_heatmap(mtx = mtx, tracks = tracks))
  expect_replays(plot_heatmap(mtx = mtx, grid = 'greens', rescale = 'cols', clust = c("average", NA)))
  expect_replays(plot_heatmap(mtx = mtx, trees = FALSE, label = FALSE))

  # A matrix from rbiom carries its own command.
  m <- taxa_matrix(rare5, rank = "Phylum", taxa = 3)
  p <- plot_heatmap(m, grid = "reds")
  expect_true(grepl("mtx <- taxa_matrix(biom", p$code, fixed = TRUE))
  expect_true(grepl("rarefy(biom", p$code, fixed = TRUE))
  expect_replays(p)
})
