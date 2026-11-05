
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


test_that("taxa_heatmap", {

  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum'))

  skip_on_cran()

  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum', legend = 'bottom', label = 'rows', trees = 'both'))
  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum', label = 'cols', trees = 'none'))
  expect_silent(taxa_heatmap(hmp5, rank = c('Phylum', 'Class'), tracks = c('Sex', 'Age')))
  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum', rescale = 'rows', title = FALSE))
  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum', tracks = list('Sex' = list(colors = c(m = "blue", f = "red")))))

  # Single-taxon and single-sample matrices.
  expect_silent(taxa_heatmap(hmp5, rank = 'Phylum', taxa = 1))
  expect_silent(taxa_heatmap(hmp5[1], rank = 'Phylum', taxa = 3))
  expect_silent(taxa_heatmap(hmp5[1], rank = 'Phylum', taxa = 1, tracks = "Sex"))

  # Tracks are per sample, so they always go on top.
  p <- taxa_heatmap(hmp5, rank = 'Phylum', tracks = list('Sex' = list(side = "left")))
  expect_true(grepl("top_track_1", p$code, fixed = TRUE))

  p <- taxa_heatmap(rare5, rank = 'Phylum', taxa = 3)
  expect_identical(as.character(p$cmd), 'taxa_heatmap(biom, "Phylum", 3)')
  expect_true(grepl("Rarefied Reads", p$code, fixed = TRUE))
  expect_true(grepl("mtx <- taxa_matrix(biom", p$code, fixed = TRUE))
  expect_true(grepl("Phylum Abundance", p$code, fixed = TRUE))
})


test_that("taxa_heatmap code replays", {

  skip_on_cran()

  b <- rare50[1:12]
  expect_replays(taxa_heatmap(b, rank = "Phylum", taxa = 5, tracks = "Sex"))
  expect_replays(taxa_heatmap(b, rank = "Phylum", taxa = 1))
  expect_replays(taxa_heatmap(b[1], rank = "Genus", taxa = 4))
  expect_replays(taxa_heatmap(b, rank = c("Phylum", "Class"), taxa = 3, tracks = c("Sex", "Age")))
  expect_replays(taxa_heatmap(b, rank = "Genus", taxa = 4, rescale = "cols", clust = c("average", NA), label_size = c(8, 6)))
})
